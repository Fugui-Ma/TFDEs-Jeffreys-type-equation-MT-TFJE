# -*- coding: utf-8 -*-
"""
%**************************************************************************
%**    Program： Robust P-dictionary multi-term CTRW; Monte Carlo MSD on [1e-6, 1e6].
%
%Explainations:
%                this expriment results showing in  Figure 1
%***Author：
%           F.G. MA
%*****Finishi Data：
%           2026-10-17
%**************************************************************************
Dependencies: numpy, torch, mpmath, matplotlib.
For admissible parameters, change ONLY P; rate and CDF grids are adaptive.
20,000 trajectories are simulated independently for each of the two laws.
Gaussian jumps are aggregated exactly in law between observation times.
"""
from pathlib import Path
from math import gamma as Gamma, pi
import csv, hashlib, json, math, os, time, warnings
import numpy as np
import torch
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

# ================= CHANGE ONLY THIS DICTIONARY =================
P = dict(alpha=0.25, beta=0.55, gamma=0.75, a=1., b=.5,
         alpha_k=(.10, .20), a_k=(.35, .20),
         beta_j=(.01, .05), b_j=(.25, .15))

# Settings unrelated to model parameters
N_PATHS, DIM, D, SEED = 20000, 1, 1., 2026
T_MIN, T_MAX, N_TIMES = 1e-6, 1e6, 37
N_THEORY, BATCH, BLOCK = 181, 256, 2048
MAX_JUMPS_PER_PATH, MAX_TOTAL_JUMPS = 2_000_000, 3_000_000_000
OUTPUT = Path(__file__).resolve().parent / 'results' if '__file__' in globals() else Path.cwd() / 'results'
FLOAT = torch.float64


class ModelError(ValueError): pass
class BudgetExceeded(RuntimeError): pass


def normalize(p):
    """Normalize input types, variable K/J, and coefficients; no hard-coded term counts."""
    req = ('alpha', 'beta', 'gamma', 'a', 'b')
    optional = ('alpha_k', 'a_k', 'beta_j', 'b_j')
    bad = set(p) - set(req) - set(optional)
    if bad: raise ModelError(f'Unknown keys in P: {sorted(bad)}')
    try:
        q = {k: float(p[k]) for k in req}
        q.update({k: tuple(float(x) for x in p.get(k, ())) for k in optional})
    except (KeyError, ValueError, TypeError) as e:
        raise ModelError('P must contain numeric alpha,beta,gamma,a,b and numeric term sequences.') from e
    if not all(math.isfinite(q[k]) for k in req) or not all(math.isfinite(x) for k in optional for x in q[k]):
        raise ModelError('P contains NaN or infinite parameters.')
    if not (0 < q['alpha'] <= 1 and 0 < q['beta'] <= 1 and 0 < q['gamma'] <= 1):
        raise ModelError('Fractional orders alpha,beta,gamma must belong to (0,1].')
    if q['a'] < 0 or q['b'] < 0: raise ModelError('Coefficients a and b must be nonnegative.')
    if len(q['alpha_k']) != len(q['a_k']) or len(q['beta_j']) != len(q['b_j']):
        raise ModelError('Lengths of alpha_k/a_k or beta_j/b_j do not match.')
    for x,v in zip(q['alpha_k'],q['a_k']):
        if not (0 < x < q['alpha'] and v >= 0): raise ModelError('Each alpha_k must lie in (0,alpha), with a_k >= 0.')
    for x,v in zip(q['beta_j'],q['b_j']):
        if not (0 < x < q['beta'] and v >= 0): raise ModelError('Each beta_j must lie in (0,beta), with b_j >= 0.')
    def poly(main_order, main_coeff, orders, coeffs):
        d = {0.: 1.}
        for exp, val in [(main_order,main_coeff), *zip(orders,coeffs)]:
            if val: d[exp] = d.get(exp,0.) + val
        return sorted(d.items())
    ap = poly(q['alpha'],q['a'],q['alpha_k'],q['a_k'])
    bp = poly(q['beta'],q['b'],q['beta_j'],q['b_j'])
    pure = len(ap)==len(bp) and all(abs(x-y)<1e-12 and abs(v-w)<1e-12 for (x,v),(y,w) in zip(ap,bp))
    linear_compound = (not q['alpha_k'] and not q['beta_j'] and
                       abs(q['alpha']-1)<1e-12 and abs(q['beta']-1)<1e-12 and
                       abs(q['gamma']-1)<1e-12 and 0 <= q['a'] <= q['b'] and q['b']>0)
    if not (pure or linear_compound or (q['beta']<=q['gamma'] and q['alpha']+q['gamma']<=1+1e-12)):
        raise ModelError('P is outside the supported admissible domain. General multi-term models require '
                         'beta <= gamma and alpha+gamma <= 1. Exact pure-power and linear compound-Poisson '
                         'exceptions are also supported; arbitrary values need not define waiting-time PDFs.')
    if N_PATHS<1 or DIM<1 or D<=0 or not (0<T_MIN<T_MAX) or N_TIMES<2:
        raise ValueError('Invalid simulation settings.')
    q['_kind'] = 'pure_power' if pure else 'linear_compound' if linear_compound else 'multiterm'
    q['_alpha_max'],q['_alpha_lead'] = ap[-1]
    q['_beta_max'],q['_beta_lead'] = bp[-1]
    delta = q['gamma']+q['_alpha_max']-q['_beta_max']
    if not (0<delta<=1+1e-10):
        raise ModelError(f'Short-time exponent delta={delta:.6g} is not in (0,1].')
    q['_delta'] = delta
    return q


def eta(z,p):
    A = 1 + p['a']*z**p['alpha'] + sum(v*z**x for x,v in zip(p['alpha_k'],p['a_k']))
    B = 1 + p['b']*z**p['beta'] + sum(v*z**x for x,v in zip(p['beta_j'],p['b_j']))
    return z**p['gamma']*A/B


def mp_model(p):
    return {k: tuple(mp.mpf(str(x)) for x in v) if isinstance(v,tuple) else mp.mpf(str(v))
            for k,v in p.items() if not k.startswith('_')}


def uniform(shape): return torch.rand(shape,dtype=FLOAT).clamp_(1e-15,1-1e-15)


def positive_stable(alpha, shape):
    """Exact positive alpha-stable: E[exp(-s S)]=exp(-s**alpha)."""
    if alpha >= 1-1e-12: return torch.ones(shape,dtype=FLOAT)
    angle = pi*uniform(shape)
    expo = -torch.log(uniform(shape))
    a = alpha
    log_s = (torch.log(torch.sin(a*angle).clamp_min(1e-300))
             - torch.log(torch.sin(angle).clamp_min(1e-300))/a
             + ((1-a)/a)*(torch.log(torch.sin((1-a)*angle).clamp_min(1e-300)) - torch.log(expo)))
    return torch.exp(log_s.clamp(max=710))


class PurePowerSampler:
    def __init__(self,p,law): self.g,self.law=p['gamma'],law
    def draw(self,shape):
        s=positive_stable(self.g,shape)
        if self.law=='rational': s=s*(-torch.log(uniform(shape)))**(1/self.g)
        return s


class LinearCompoundSampler:
    """Exact for eta(z)=z(1+a*z)/(1+b*z), 0<=a<=b."""
    def __init__(self,p,law):
        self.drift=p['a']/p['b']; self.intensity=(1-self.drift)/p['b']
        self.scale=p['b']; self.law=law
    def draw(self,shape):
        # NumPy generator is seeded by torch initial_seed; draws are reproducible for each law.
        rng=np.random.default_rng(int(torch.randint(0,2**31-1,(1,)).item()))
        n=int(np.prod(shape))
        u=np.ones(n) if self.law=='exponential' else rng.exponential(1,size=n)
        k=rng.poisson(self.intensity*u)
        jumps=np.zeros(n,dtype=float); active=k>0
        if np.any(active): jumps[active]=rng.gamma(shape=k[active],scale=self.scale)
        return torch.from_numpy((self.drift*u+jumps).reshape(shape).copy())


class RationalMixSampler:
    """Adaptive nonnegative Stieltjes spectral mixture of exponentials."""
    def __init__(self,p):
        self.p=p
        for decades in (32,48,64,88,120,160,210,260):
            n=int(450*2*decades)+1
            u=torch.linspace(-decades*math.log(10),decades*math.log(10),n,dtype=FLOAT)
            r=torch.exp(u)
            z=torch.polar(r,torch.full_like(r,math.pi))
            w=-(1/(1+eta(z,p))).imag/math.pi
            if not torch.isfinite(w).all() or w.min().item() < -2e-10:
                raise ModelError('P fails the nonnegative Stieltjes spectral density check; '
                                 'the proposed waiting-time distribution may be invalid.')
            w=w.clamp_min(0)
            area=(w[:-1]+w[1:])*(u[1:]-u[:-1])/2
            mass=float(area.sum())
            if abs(mass-1) > 3e-5: continue
            good=True
            for s in (1/T_MAX,1.,1/T_MIN):
                check=float(torch.trapezoid(w*r/(r+s),u))
                theoretical=1/(1+float(eta(s,p)))
                if abs(check-theoretical)>5e-5: good=False; break
            if good:
                self.log_rates=u
                self.cdf=torch.cat([torch.zeros(1,dtype=FLOAT),area.cumsum(0)/mass]);self.cdf[-1]=1.
                print(f'  Rational mixing range: +/-{decades} decades, mass={mass:.10f}',flush=True)
                return
        raise ModelError('This P requires a wider spectral sampling range than float64 can resolve '
                         '(possible extremely small delta). No inaccurate waiting-time samples were used.')
    def draw(self,shape):
        u=uniform(shape).reshape(-1)
        k=torch.searchsorted(self.cdf,u,right=True).clamp(1,len(self.cdf)-1)
        f=(u-self.cdf[k-1])/(self.cdf[k]-self.cdf[k-1]).clamp_min(1e-300)
        rate=torch.exp(self.log_rates[k-1]+f*(self.log_rates[k]-self.log_rates[k-1]))
        return (-torch.log(uniform(rate.shape))/rate).reshape(shape)


class ExponentialCDFSampler:
    """Parameter-specific, adaptively verified CDF for psi_hat=exp(-eta)."""
    def __init__(self,p):
        OUTPUT.mkdir(parents=True,exist_ok=True)
        self.p=p
        key=hashlib.sha256(json.dumps({'p':{k:v for k,v in p.items() if not k.startswith('_')},
                                       'time':[T_MIN,T_MAX],'version':'adaptive-cdf-4'},sort_keys=True).encode()).hexdigest()
        cache=OUTPUT/f'exponential_cdf_{key[:16]}.npz'
        if cache.exists():
            with np.load(cache,allow_pickle=False) as obj:
                self.F=torch.tensor(obj['F'].copy(),dtype=FLOAT)
                self.logt=torch.tensor(obj['logt'].copy(),dtype=FLOAT)
            self.verify(p)
            print('  Exponential CDF cache matches current P',flush=True)
            return
        self.pm=mp_model(p)
        self.calls=0
        grid=np.geomspace(T_MIN/1e4,T_MAX,110)
        F=self.evaluate(grid)
        for _ in range(12):
            if F[0]<1e-8: break
            g=grid[0]/100
            grid=np.r_[g,grid];F=np.r_[self.evaluate([g]),F]
        if F[0]>1e-8:
            raise ModelError('Exponential waiting-time distribution has nonnegligible mass below '
                             'adaptive CDF lower range; cannot safely interpolate.')
        # Only refine parts of the CDF that carry probability.
        for iteration in range(12):
            difference=np.diff(F)
            if np.min(difference)<-3e-7:
                raise ModelError('Computed CDF is nonmonotone for P; increase precision or use a different law.')
            midpoint=np.sqrt(grid[:-1]*grid[1:])
            active=np.where((difference>0.007) | ((difference>1e-5) & (np.log(grid[1:]/grid[:-1])>.17)))[0]
            if not active.size: break
            if len(grid)+len(active)>1600:
                raise ModelError('Adaptive exponential waiting-time CDF exceeds 1600 nodes; distribution too sharp.')
            newt=midpoint[active]
            newF=self.evaluate(newt)
            grid=np.r_[grid,newt]; F=np.r_[F,newF]
            order=np.argsort(grid);grid,F=grid[order],F[order]
        if np.min(np.diff(F))<-3e-7:
            raise ModelError('Exponential CDF precision insufficient for current P.')
        if np.max(np.diff(F)) > 0.008:
            raise ModelError('Adaptive exponential CDF did not reach its target resolution; no unverified samples were drawn.')
        F=np.maximum.accumulate(F)
        keep=np.r_[True,np.diff(F)>1e-13]
        grid,F=grid[keep],F[keep]
        if len(F)<2: raise ModelError('Exponential CDF contains no resolved mass on the time interval.')
        self.logt=torch.tensor(np.log(grid),dtype=FLOAT)
        self.F=torch.tensor(F,dtype=FLOAT)
        self.verify(p)
        np.savez_compressed(cache,logt=self.logt.numpy(),F=self.F.numpy())
        print(f'  Exponential adaptive CDF: {len(F)} nodes, {self.calls} inverse transforms, F(tmax)={F[-1]:.7f}',flush=True)
    def evaluate(self,t):
        values=[]
        for x in t:
            # Rigorous one-sided Chernoff bound: F(t)<=exp(s*t-eta(s)).
            # Prevents invalid Talbot inversions of exponentially small CDF values.
            xx=float(x)
            z=(np.array([1,2,4,8,16,32,64,128],dtype=float)/xx)
            eta_val=np.asarray(eta(z,self.p),dtype=float)
            score=np.max(eta_val-z*xx)
            if score>26:  # F(t)<exp(-26)<5.2e-12; negligible for MC.
                values.append(0.);self.calls+=1;continue
            tt=mp.mpf(str(xx))
            with mp.workdps(65):
                val=mp.invertlaplace(lambda z: mp.exp(-eta(z,self.pm))/z,tt,method='talbot',degree=90)
            # Fixed Talbot can fail very near the super-exponentially small tail.
            if not mp.isfinite(val) or val < -mp.mpf('1e-12') or val > 1+mp.mpf('1e-12') or (val<mp.mpf('1e-20')):
                with mp.workdps(75):
                    val=mp.invertlaplace(lambda z: mp.exp(-eta(z,self.pm))/z,tt,method='dehoog',degree=65)
            if not mp.isfinite(val) or val < -mp.mpf('1e-9') or val > 1+mp.mpf('1e-9'):
                raise ModelError(f'Exponential waiting CDF cannot be inverted reliably at t={xx:.5g}; '
                                 'the distribution may have an atom or numerical conditioning may be poor.')
            values.append(float(min(1,max(0,val))))
            self.calls+=1
        return np.asarray(values)
    def verify(self,p):
        F=self.F.numpy();logt=self.logt.numpy()
        if len(F)<2 or np.any(np.diff(F)<=0) or F[0]>1e-8:
            raise ModelError('Exponential waiting-time CDF is not strictly increasing or has missing lower mass.')
        # Independent transform reconstruction using Simpson quadrature on log-time interpolation.
        left=np.exp(logt[:-1]);right=np.exp(logt[1:]);mid=np.sqrt(left*right)
        for s in (1/T_MAX,1e-3,1.,1e3,1/T_MIN):
            expected=math.exp(-float(eta(s,p)))
            if expected < 1e-6: continue
            # Laplace-Stieltjes transform approximation on the CDF mesh.
            approximate=np.sum(np.diff(F)*(np.exp(-s*left)+4*np.exp(-s*mid)+np.exp(-s*right))/6)
            # Upper tail beyond T_MAX contributes at most 1-F(T_MAX).
            error_bound=(1-F[-1])*np.exp(-s*right[-1])
            if abs(approximate-expected)>0.008+0.006*expected+error_bound:
                raise ModelError(f'Exponential CDF transform inconsistent with P at s={s:g}: '
                                 f'{approximate:.6g} vs {expected:.6g}.')
    def draw(self,shape):
        u=uniform(shape).reshape(-1)
        k=torch.searchsorted(self.F,u,right=True).clamp(1,len(self.F)-1)
        v=(u-self.F[k-1])/(self.F[k]-self.F[k-1]).clamp_min(1e-300)
        out=torch.exp(self.logt[k-1]+v*(self.logt[k]-self.logt[k-1]))
        out[u>=self.F[-1]]=float('inf')   # exact censoring through t<=T_MAX
        return out.reshape(shape)


def sampler(p,law):
    if p['_kind']=='pure_power': return PurePowerSampler(p,law)
    if p['_kind']=='linear_compound': return LinearCompoundSampler(p,law)
    return RationalMixSampler(p) if law=='rational' else ExponentialCDFSampler(p)


def validate_sampler(p,law,obj):
    torch.manual_seed(SEED+170)
    waiting=obj.draw((50000,))
    if torch.isnan(waiting).any() or (waiting<0).any(): raise ModelError(f'{law} sampler emitted invalid waiting times.')
    for s in (1/T_MAX,1.,1/T_MIN):
        exact=1/(1+float(eta(s,p))) if law=='rational' else math.exp(-float(eta(s,p)))
        if exact<.02: continue
        estimate=float(torch.exp(-s*waiting).mean())
        if abs(estimate-exact)>max(.018,.04*exact):
            raise ModelError(f'{law} waiting-time Monte Carlo Laplace check failed at s={s:g}: '
                             f'{estimate:.5f} vs {exact:.5f}.')


def integer_power_monte_carlo(obs,npaths,law,seed,save_paths=True):
    """Exact renewal-count sampler for eta(s)=s: Poisson or deterministic."""
    torch.manual_seed(seed)
    t=torch.as_tensor(obs,dtype=FLOAT)
    difft=torch.diff(torch.cat([torch.zeros(1,dtype=FLOAT),t]))
    ntimes=len(obs);s1=torch.zeros(ntimes,dtype=FLOAT);s2=torch.zeros_like(s1)
    nonzero=torch.zeros(ntimes,dtype=torch.int64);sumN=torch.zeros(ntimes,dtype=FLOAT)
    pos_all=[];n_all=[]
    for start in range(0,npaths,BATCH):
        b=min(BATCH,npaths-start)
        if law=='rational':
            inc=torch.poisson(difft.expand(b,-1).contiguous()).to(torch.int64)
            counts=inc.cumsum(1)
        else:
            counts=torch.floor(t*(1+1e-14)).to(torch.int64).expand(b,-1).clone()
            inc=torch.diff(torch.cat([torch.zeros((b,1),dtype=torch.int64),counts],dim=1),dim=1)
        steps=torch.randn((b,ntimes,DIM),dtype=FLOAT)*torch.sqrt(2*D*inc.to(FLOAT))[:,:,None]
        pos=steps.cumsum(1)
        sq=(pos*pos).sum(2)
        s1+=sq.sum(0);s2+=(sq*sq).sum(0);nonzero+=(counts>0).sum(0);sumN+=counts.sum(0)
        if save_paths: pos_all.append(pos.numpy().copy());n_all.append(counts.numpy().copy())
    mse=s1.numpy()/npaths
    se=np.sqrt(torch.clamp((s2-s1*s1/npaths)/(npaths-1),min=0).numpy()/npaths) if npaths>1 else np.full(ntimes,np.nan)
    print(f'  {law}: exact integer-power renewal counts, {npaths:,} paths',flush=True)
    return dict(msd=mse,se=se,nonzero=nonzero.numpy(),meanN=sumN.numpy()/npaths,
                positions=np.concatenate(pos_all) if save_paths else None,
                counts=np.concatenate(n_all) if save_paths else None)


def monte_carlo(p,law,obs,npaths=N_PATHS,seed=SEED,save_paths=True,obj=None):
    if p['_kind']=='pure_power' and abs(p['gamma']-1)<1e-12:
        return integer_power_monte_carlo(obs,npaths,law,seed,save_paths)
    if obj is None: obj=sampler(p,law)
    validate_sampler(p,law,obj)
    torch.manual_seed(seed)
    t=torch.as_tensor(obs,dtype=FLOAT)
    m=t.numel();s1=torch.zeros(m,dtype=FLOAT);s2=torch.zeros_like(s1)
    nonzero=torch.zeros(m,dtype=torch.int64);sumN=torch.zeros(m,dtype=FLOAT)
    positions_all=[];counts_all=[];done=0;event_total=0;tic=time.monotonic()
    while done<npaths:
        B=min(BATCH,npaths-done)
        clock=torch.zeros(B,dtype=FLOAT);ncount=torch.zeros(B,dtype=torch.int64)
        nt=torch.zeros((B,m),dtype=torch.int64)
        decided=torch.zeros((B,m),dtype=torch.bool)
        active=torch.arange(B)
        while active.numel():
            if int(ncount[active].max())+BLOCK>MAX_JUMPS_PER_PATH:
                raise BudgetExceeded(f'{law}: some paths require over {MAX_JUMPS_PER_PATH:,} jumps by t={T_MAX:g} '
                                     'for this P. Increase budget or reduce horizon; no paths were truncated.')
            w=obj.draw((active.numel(),BLOCK))
            end=clock[active,None]+w.cumsum(dim=1)
            crossed=(~decided[active])&(t[None,:]<=end[:,-1,None])
            if crossed.any():
                k=torch.searchsorted(end.contiguous(),t.expand(len(active),-1).contiguous(),right=True)
                vals=ncount[active,None]+k
                nt[active]=torch.where(crossed,vals,nt[active]);decided[active]|=crossed
            clock[active]=end[:,-1];ncount[active]+=BLOCK
            active=active[clock[active]<t[-1]]
        if not decided.all() or not (nt[:,1:]>=nt[:,:-1]).all(): raise RuntimeError('Renewal recording logic failed.')
        diffN=torch.diff(torch.cat((torch.zeros((B,1),dtype=torch.int64),nt),dim=1),dim=1)
        dx=torch.randn((B,m,DIM),dtype=FLOAT)*torch.sqrt(2*D*diffN.to(FLOAT))[:,:,None]
        pos=dx.cumsum(dim=1)
        r2=(pos*pos).sum(dim=2)
        s1+=r2.sum(0);s2+=(r2*r2).sum(0);nonzero+=(nt>0).sum(0);sumN+=nt.sum(0)
        if save_paths: positions_all.append(pos.numpy().copy());counts_all.append(nt.numpy().copy())
        done+=B;event_total+=int(nt[:,-1].sum())
        if done%(BATCH*10)==0 or done==npaths:
            print(f'  {law}: {done:,}/{npaths:,}; avg renewals={event_total/done:,.0f}; '
                  f'elapsed={time.monotonic()-tic:.1f}s',flush=True)
    mean=s1.numpy()/npaths
    se=np.sqrt(torch.clamp((s2-s1*s1/npaths)/(npaths-1),min=0).numpy()/npaths) if npaths>1 else np.full(m,np.nan)
    return dict(msd=mean,se=se,nonzero=nonzero.numpy(),meanN=sumN.numpy()/npaths,
                positions=np.concatenate(positions_all) if save_paths else None,
                counts=np.concatenate(counts_all) if save_paths else None)


def fixed_talbot(F,t,N=32):
    x=torch.tensor(np.asarray(t),dtype=FLOAT);theta=pi*torch.arange(1,N,dtype=FLOAT)/N
    r=.4*N;ct=1/torch.tan(theta)
    z0=(r/x).to(torch.complex128)
    z=(r/x[:,None])*theta[None,:]*(ct[None,:]+1j)
    w=1+1j*theta*(1+ct**2)-1j*ct
    return ((.4/x)*(.5*math.exp(r)*F(z0).real+(torch.exp(x[:,None]*z)*w[None,:]*F(z)).real.sum(1))).numpy()


def theoretical_msd(p,t,law,reference=None):
    def fun(z):
        h=eta(z,p)
        if law=='rational': q=1/h
        else:
            q=torch.empty_like(h);positive=h.real>=0
            q[positive]=torch.exp(-h[positive])/(-torch.expm1(-h[positive]))
            q[~positive]=1/torch.expm1(h[~positive])
        return 2*DIM*D*q/z
    x=np.asarray(t,dtype=float)
    if p['_kind']=='pure_power' and abs(p['gamma']-1)<1e-12:
        r=2*DIM*D*x if law=='rational' else 2*DIM*D*np.floor(x)
        zero=(r==0)
        r=r.astype(float);r[zero]=np.nan
        return r,zero
    lo,hi=fixed_talbot(fun,x,24),fixed_talbot(fun,x,32)
    error=np.abs(hi-lo)/np.maximum(abs(hi),1e-290)
    invalid=(~np.isfinite(hi))|(hi<=0)|(~np.isfinite(error))|(error>3e-7)
    censored=np.zeros(len(x),dtype=bool)
    if law=='exponential':
        if reference is None: raise ValueError('Exponential MSD needs rational reference for safe short-time bounds.')
        # Chernoff bound optimized over several positive Laplace arguments s:
        # E[N(t)] <= exp(s*t)/(exp(eta(s))-1), valid for all renewal processes.
        # Do not invert values that are provably negligible relative to rational MSD.
        factors=np.array([1.,2.,4.,8.,16.,32.,64.,128.,256.,512.])
        s=factors[:,None]/x[None,:]
        h=np.asarray(eta(s,p),dtype=float)
        with np.errstate(over='ignore',under='ignore',divide='ignore',invalid='ignore'):
            logden=np.empty_like(h)
            large=h>1
            logden[large]=h[large]+np.log1p(-np.exp(-h[large]))
            logden[~large]=np.log(np.expm1(h[~large]))
        logupper=np.min(s*x[None,:]-logden,axis=0)+math.log(2*DIM*D)
        censored=(logupper-np.log(reference))<math.log(1e-9)
        hi[censored]=np.nan;invalid[censored]=False
        invalid |= (x<1e-5)&(~censored)
    pp=mp_model(p)
    for j in np.flatnonzero(invalid):
        val=None
        for digits,degree,method in ((60,90,'talbot'),(90,75,'dehoog')):
            with mp.workdps(digits):
                def fm(z):
                    v=eta(z,pp);q=1/v if law=='rational' else 1/mp.expm1(v)
                    return 2*DIM*mp.mpf(str(D))*q/z
                result=mp.invertlaplace(fm,mp.mpf(str(float(x[j]))),method=method,degree=degree)
            if mp.isfinite(result) and result>0 and float(result)>0:
                val=float(result);break
        if val is None or val==0: raise ArithmeticError(f'Theoretical {law} MSD could not be inverted at t={x[j]:g}.')
        hi[j]=val
    if not np.all(np.isfinite(hi[~censored])&(hi[~censored]>0)):
        raise ArithmeticError(f'Nonpositive theoretical MSD ({law}).')
    return hi,censored




# ============================================================
# Plot settings
# ============================================================
PLOT = dict(
    figsize=(8.6, 5.5),
    axis_font=24, tick_font=18, legend_font=10, asym_font=20,
    curve_lw=3.2, guide_lw=2.8, E_lw=2.8,
    marker_size=8, marker_edge=1.0, legend_cols=4,
    line_styles=(
    '-',                      # Rational MSD
    '--',                     # Exponential MSD
    ':',                      # Short-time asymptotic
    '-.',                     # Long-time asymptotic
    (0, (8, 2, 1, 2, 1, 2))  # E(t): custom dash pattern
     ),
    marker_shapes=('o', '^'),
    min_nonzero=20, max_rel_se=.45,
    errorbars=False, dpi=300,
    colors=('#1f77b4', '#ff7f0e', '#2ca02c',
            '#d62728', '#7b5bbd')
)


# ============================================================
# Log-log MSD plot with two-column legend
# ============================================================
def plot_result(p, t, rat, exp, cens, obs, rm, em,
                npaths=N_PATHS, style=None):

    st = {**PLOT, **(style or {})}
    blue, orange, green, red, purple = st['colors']
    p = p if '_delta' in p else normalize(p)

    t = np.asarray(t, float)
    obs = np.asarray(obs, float)
    rat = np.asarray(rat, float)
    exp = np.asarray(exp, float)
    cens = np.asarray(cens, bool)

    if not (len(t) == len(rat) == len(exp) == len(cens)
            and len(obs) == len(rm['msd']) == len(em['msd'])):
        raise ValueError('Time and MSD array lengths differ')

    if not (np.all(np.diff(t) > 0)
            and np.all(np.diff(obs) > 0)
            and np.all(rat > 0)):
        raise ValueError('Invalid time grid or rational MSD')

    # Asymptotic laws
    delta, g = p['_delta'], p['gamma']
    short = (2*DIM*D*p['_beta_lead']/p['_alpha_lead']
             /Gamma(1+delta)*t**delta)
    long = 2*DIM*D/Gamma(1+g)*t**g

    mshort = t <= t[0]*(t[-1]/t[0])**.32
    mlong = t >= t[0]*(t[-1]/t[0])**.78

    # Relative discrepancy
    E = np.full_like(rat, np.nan)
    E[cens] = 1.0
    good = (~cens) & np.isfinite(exp)
    E[good] = np.abs(rat[good]-exp[good])/rat[good]

    # Monte Carlo filtering
    def keep(data):
        y = np.asarray(data['msd'], float)
        mask = np.isfinite(y) & (y > 0)
        if 'nonzero' in data:
            mask &= np.asarray(data['nonzero']) >= st['min_nonzero']
        if 'se' in data:
            mask &= np.asarray(data['se'])/np.maximum(y, 1e-300) <= st['max_rel_se']
        return mask

    kr, ke = keep(rm), keep(em)

    # Automatic y-axis range
    vals = np.r_[rat, np.asarray(rm['msd'])[kr],
                 np.asarray(em['msd'])[ke]]
    ve = exp[(~cens) & np.isfinite(exp) & (exp > rat.min()*1e-4)]
    if len(ve):
        vals = np.r_[vals, ve]

    low = max(1e-14, vals.min()/4, rat.min()*1e-4)
    high = vals.max()*1.8

    with plt.rc_context({
        'font.family': 'serif',
        'mathtext.fontset': 'stix',
        'font.size': st['tick_font'],
        'savefig.dpi': st['dpi']
    }):

        fig, ax = plt.subplots(figsize=st['figsize'])
        ax2 = ax.twinx()

        # Frame
        for side in ('left', 'bottom', 'top'):
            ax.spines[side].set_visible(True)
            ax.spines[side].set_linewidth(1.1)
            ax2.spines[side].set_visible(False)

        ax.spines['right'].set_visible(False)
        ax2.spines['right'].set_visible(True)
        ax2.spines['right'].set_color(purple)

        # Theoretical MSDs
        r, = ax.loglog(t, rat, '-', color=blue,
                       lw=st['curve_lw'], label='$M_{\mathrm{rat}}$')
        e, = ax.loglog(t, exp, '--', color=orange,
                       lw=st['curve_lw'], label='M_{\mathrm{exp}}')

        # Asymptotic lines
        s, = ax.loglog(t[mshort], short[mshort], ':',
                       color='black', lw=st['guide_lw'],
                       label=rf'Short $\sim t^{{{delta:g}}}$')
        l, = ax.loglog(t[mlong], long[mlong], '-.',
                       color=red, lw=st['guide_lw'],
                       label=rf'Long $\sim t^{{{g:g}}}$')

        # Relative discrepancy
        ed, = ax2.loglog(t, np.where(E > 0, E, np.nan),
                         '-.', color=purple, lw=st['E_lw'],
                         label=r'$E(t)=|M_{\mathrm{rat}}-M_{\mathrm{exp}}|/M_{\mathrm{rat}}$')

        # Monte Carlo points
        pr, = ax.loglog(
            obs[kr], np.asarray(rm['msd'])[kr],
            ls='none', marker='o', ms=st['marker_size'],
            mfc='white', mec=blue, mew=st['marker_edge'],
            label='Rational MC'
        )

        pe, = ax.loglog(
            obs[ke], np.asarray(em['msd'])[ke],
            ls='none', marker='^', ms=st['marker_size'],
            mfc='white', mec=orange, mew=st['marker_edge'],
            label='Exponential MC'
        )

        # Optional error bars
        if st['errorbars']:
            for data, mask, color in ((rm, kr, blue), (em, ke, orange)):
                if 'se' in data and np.any(mask):
                    y = np.asarray(data['msd'])[mask]
                    se = np.asarray(data['se'])[mask]
                    ax.errorbar(
                        obs[mask], y,
                        yerr=np.minimum(se, .9*y),
                        fmt='none', color=color,
                        alpha=.55, capsize=2, lw=.8
                    )

        # Axis limits
        ax.set(xlim=(t[0], t[-1]), ylim=(low, high))
        pos = E[np.isfinite(E) & (E > 0)]
        ax2.set_ylim(
            max(1e-12, min(20e-6, pos.min()/2)) if len(pos) else 20e-6,
            max(20, 20*pos.max()) if len(pos) else 20
        )

        # Axis labels
        ax.set_xlabel(r'Time $t$ (log scale)',
                      fontsize=st['axis_font'], labelpad=8)
        ax.set_ylabel(r'MSD $M(t)$ (log scale)',
                      fontsize=st['axis_font'], labelpad=8)
        ax2.set_ylabel(r'$E(t)$ (log scale)',
                       fontsize=st['axis_font'],
                       color=purple, labelpad=8)

        # Ticks
        for a in (ax.xaxis, ax.yaxis, ax2.yaxis):
            a.set_minor_locator(
                LogLocator(base=10, subs=np.arange(2, 10), numticks=200)
            )
            a.set_minor_formatter(NullFormatter())

        ax.tick_params(axis='x', which='both', top=True,
                       direction='in', labelsize=st['tick_font'])
        ax.tick_params(axis='y', which='both', left=True,
                       direction='in', labelsize=st['tick_font'])
        ax2.tick_params(axis='y', which='both', right=True,
                        direction='in', colors=purple,
                        labelsize=st['tick_font'])
        ax.grid(which='major', alpha=.15)

        # Asymptotic exponent labels
        def annotate_power(xx, yy, exponent, color, early):
            ii = np.flatnonzero(
                (yy >= low*6) & (yy <= high/6) & np.isfinite(yy)
            )

            label = rf'$\sim t^{{{exponent:g}}}$'

            if len(ii):
                j = ii[int((.40 if early else .40)*(len(ii)-1))]
                ax.annotate(
                    label, (xx[j], yy[j]),
                    xytext=(8, -20) if early else (-12, 15),
                    textcoords='offset points',
                    fontsize=st['asym_font'], color=color,
                    ha='left' if early else 'right',
                    va='top' if early else 'bottom'
                )
            else:
                ax.text(
                    .12 if early else .70,
                    .24 if early else .80,
                    label, transform=ax.transAxes,
                    fontsize=st['asym_font'], color=color,
                    ha='left',
                    va='center',
                    zorder=20
                )

        annotate_power(t[mshort], short[mshort], delta, 'black', True)
        annotate_power(t[mlong], long[mlong], g, red, False)

        # ===================================================
        # Legend: INSIDE, lower right, TWO COLUMNS
        # ===================================================
        hs = [r, pr, s, ed, e, pe, l]
        labels = [
            '$M_{\mathrm{rat}}$', 'Rational MC',
            rf'Short time $\sim t^{{{delta:g}}}$', r'$E(t)=|M_{\mathrm{rat}}-M_{\mathrm{exp}}|/M_{\mathrm{rat}}$',
            '$M_{\mathrm{exp}}$', 'Exponential MC',
            rf'Long time $\sim t^{{{g:g}}}$'
        ]

        for h, text in zip(hs, labels):
            h.set_label(text)

        ax.set_zorder(ax2.get_zorder() + 1)
        ax.patch.set_visible(False)

        ax.legend(
            handles=hs,
            loc='lower right',
            bbox_to_anchor=(0.975, 0.045),
            bbox_transform=ax.transAxes,
            ncol=2,
            fontsize=9,
            frameon=True,
            framealpha=.96,
            facecolor='white',
            edgecolor='0.75',
            columnspacing=.9,
            labelspacing=.35,
            handlelength=1.8,
            handletextpad=.45,
            borderpad=0.55,
            borderaxespad=0.4
        )
        
        fig.subplots_adjust(
            left=.17, right=.82, bottom=.19, top=.96
        )

        OUTPUT.mkdir(parents=True, exist_ok=True)

        for ext in ('png', 'pdf'):
            fig.savefig(
                OUTPUT/f'Figure1_msd_a.{ext}',
                bbox_inches='tight'
            )

        plt.close(fig)

    return E


def main(p=None,npaths=N_PATHS,output=None,do_theory=True,save_paths=True):
    global OUTPUT
    if npaths<1: raise ValueError('npaths must be >=1')
    p=normalize(P if p is None else p)
    directory=Path(OUTPUT if output is None else output);directory.mkdir(parents=True,exist_ok=True)
    OUTPUT=directory
    torch.set_num_threads(min(4,os.cpu_count() or 1))
    obs=np.geomspace(T_MIN,T_MAX,N_TIMES)
    start=time.monotonic()
    print(f"P classification: {p['_kind']}; gamma={p['gamma']}, short-time delta={p['_delta']:.6g}")
    theory_t=np.geomspace(T_MIN,T_MAX,N_THEORY)
    theory_r=theory_e=cens=None
    if do_theory:
        theory_r,_=theoretical_msd(p,theory_t,'rational')
        theory_e,cens=theoretical_msd(p,theory_t,'exponential',reference=theory_r)
        # A transparent preflight: do not start billions of events if resource budget is exceeded.
        mu=float(theory_r[-1]/(2*DIM*D))
        predicted=2*npaths*mu
        print(f'Estimated total waiting-time draws (both laws): {predicted:.3g}')
        if predicted>MAX_TOTAL_JUMPS and not (p['_kind']=='pure_power' and abs(p['gamma']-1)<1e-12):
            raise BudgetExceeded(f'For P the expected total number of simulated waits is {predicted:.3g}, '
                                 f'greater than budget {MAX_TOTAL_JUMPS:.3g}. Cannot guarantee 20,000 paths '
                                 'through t=1e6 at the present computational budget.')
    # Prepare BOTH waiting-time generators before starting expensive simulations.
    print('Building and validating both waiting-time samplers...',flush=True)
    sampler_rat=sampler(p,'rational')
    sampler_exp=sampler(p,'exponential')
    validate_sampler(p,'rational',sampler_rat)
    validate_sampler(p,'exponential',sampler_exp)
    rat=monte_carlo(p,'rational',obs,npaths,seed=SEED,save_paths=save_paths,obj=sampler_rat)
    exp=monte_carlo(p,'exponential',obs,npaths,seed=SEED+333,save_paths=save_paths,obj=sampler_exp)
    if do_theory: E=plot_result(p,theory_t,theory_r,theory_e,cens,obs,rat,exp,npaths=npaths)
    else: E=None
    with (directory/'mc_observations.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['t','M_rat_MC','rat_SE','rat_nonzero','rat_meanN','M_exp_MC','exp_SE','exp_nonzero','exp_meanN'])
        w.writerows(zip(obs,rat['msd'],rat['se'],rat['nonzero'],rat['meanN'],exp['msd'],exp['se'],exp['nonzero'],exp['meanN']))
    if do_theory:
        with (directory/'msd_theory.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.writer(f);w.writerow(['t','M_rat_theory','M_exp_theory','E_theory','exp_censored'])
            w.writerows(zip(theory_t,theory_r,theory_e,E,cens))
    if save_paths:
        np.savez_compressed(directory/'trajectory_observations.npz',t=obs,rat_positions=rat['positions'],
                            exp_positions=exp['positions'],rat_counts=rat['counts'],exp_counts=exp['counts'])
    (directory/'model_parameters.json').write_text(json.dumps({k:v for k,v in p.items() if not k.startswith('_')},indent=2),encoding='utf-8')
    print(f'Completed in {time.monotonic()-start:.1f}s; results in {directory}',flush=True)
    return {'times':obs,'rational':rat,'exponential':exp}


if __name__=='__main__': main()
