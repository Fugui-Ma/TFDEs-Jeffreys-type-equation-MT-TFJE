function out=J9_source_solver(M,N,times,dim,caseName,p,method)
% C9_SOURCE_SOLVER
% Core solver adapted from the uploaded Jeffreys_verify.m.
%
% It preserves the original:
%   - hyperbolic contour quadrature;
%   - Chebyshev-Gauss-Lobatto interpolation (FFT or direct DCT-I);
%   - homogeneous Dirichlet Legendre-Galerkin discretization;
%   - sparse complex-valued A\rhs solve per contour node (mldivide).
%
% Special:
%   J9_source_solver('check')  runs the mathematical consistency checks.
%
% No files are written. T.stage has EIGHT measured, nonoverlapping
% components summing to T.total from the same execution.

if nargin==1 && (ischar(M)||isstring(M))
    if strcmpi(char(M),'check')
        self_test();
        out=[];
        return
    else
        error('Unknown diagnostic action.');
    end
end

times=times(:);
assert(any(dim==[1 2]),'Dimension must be 1 or 2.');
assert(M>=3 && M==round(M),'M must be an integer >=3.');
assert(N>=2 && N==round(N),'N must be an integer >=2.');
assert(all(times>=p.t0-1e-13 & times<=p.Lambda*p.t0+1e-13), ...
    'Times must lie in [t0,Lambda*t0].');
assert(any(strcmpi(method,{'fft','direct'})), ...
    'Method must be fft or direct.');

wall=tic;

% Nine nonoverlapping timing categories:
%  1 basis/contour setup plus one-off spatial-shape or initial-data setup
%  2 evaluation of fhat (1D) or Laplace symbols/initial RHS (2D f=0)
%  3 CGL DCT-I transform
%  4 source/RHS vector assembly via Chebyshev-to-Legendre projection
%  5 per-contour-node Galerkin matrix assembly
%  6 complete A\rhs linear-system solution (MATLAB mldivide)
%  7 inverse-Laplace/CIM time accumulation
%  8 physical-space solution reconstruction
T.setup=0;
T.data=0;
T.fhat=0;
T.dct=0;
T.projection=0;
T.assembly=0;
T.linearSolve=0;
T.accum=0;
T.reconstruct=0;

% One-time spatial and contour setup.
phase=tic;
B=make_basis(M);
m=M-1;
ndof=m^dim;
if strcmpi(method,'direct')
    C=cosine_dense(M);
else
    C=[];
end

if dim==1
    SpatialMass=B.mass;
    SpatialStiff=B.stiff;
else
    SpatialMass=kron(B.mass,B.mass);
    SpatialStiff=kron(B.mass,B.stiff)+kron(B.stiff,B.mass);
end
[z,zp,tau]=hyperbola_nodes(N,p.t0,p.Lambda);
T.setup=T.setup+toc(phase);

% Prepare either the six-peak initial condition (2D) or the manufactured
% solution shape and its Laplacian (1D).
phase=tic;
if strcmpi(caseName,'sixpeak')
    assert(dim==2,'Six-peak initial value is 2D.');
    [X,Y]=ndgrid(B.nodes,B.nodes);
    initial=six_gaussians(X,Y);
    vn=[];
    lapvn=[];
else
    [v,vpp]=shape1(B.nodes,caseName);
    if dim==1
        vn=v;
        lapvn=vpp;
    else
        vn=v*v.';
        lapvn=vpp*v.'+v*vpp.';
    end
    initial=[];
end
T.data=T.data+toc(phase);

if strcmpi(caseName,'sixpeak')
    % Transform p_0 once, rather than N times.
    [b0,dctTime,projTime]=modal_load(B,initial,method,C);
    T.dct=T.dct+dctTime;
    T.projection=T.projection+projTime;
else
    b0=zeros(ndof,1);
end

coeff=zeros(numel(times),ndof);

for k=1:N
    zk=z(k);

    % Source and coefficient-symbol evaluation. For f=0, this phase only
    % evaluates the Laplace symbols and forms the initial-value RHS.
    phase=tic;
    [nz,dz]=laplace_terms(zk,p);
    if strcmpi(caseName,'sixpeak')
        rhs=(zk^(p.gamma-1)*nz)*b0;
    else
        % Manufactured forcing in the Laplace domain:
        % fhat(z) = Gamma(kappa+1)/z^(kappa+1) *
        %          [z^gamma*N(z)*v - D(z)*Laplacian(v)].
        fhat=gamma(p.kappa+1)*zk^(-p.kappa-1)* ...
            (zk^p.gamma*nz*vn-dz*lapvn);
    end
    T.fhat=T.fhat+toc(phase);

    if ~strcmpi(caseName,'sixpeak')
        % Measure transform and projection SEPARATELY.
        [rhs,dctTime,projTime]=modal_load(B,fhat,method,C);
        T.dct=T.dct+dctTime;
        T.projection=T.projection+projTime;
    end

    % One complex Galerkin matrix, new for each z_k.
    phase=tic;
    A=(zk^p.gamma*nz)*SpatialMass+dz*SpatialStiff;
    T.assembly=T.assembly+toc(phase);

    % Full MATLAB sparse linear solve A*chat = rhs.
    % mldivide (A\rhs) includes its internally chosen factorization
    % and substitution work; these are NOT timed separately.
    phase=tic;
    chat=A\rhs;
    T.linearSolve=T.linearSolve+toc(phase);

    % Evaluate and accumulate the CIM inverse-Laplace contribution.
    phase=tic;
    coeff=coeff+(tau/pi)*imag((exp(zk*times)*zp(k))*chat.');
    T.accum=T.accum+toc(phase);
end

% Reconstruct the solution from its Galerkin coefficients.
phase=tic;
xVis=linspace(-1,1,32)';
Vvis=legendre_dirichlet(xVis,M);
for k=1:numel(times)
    if dim==1
        field=Vvis*coeff(k,:).'; %#ok<NASGU>
    else
        field=Vvis*reshape(coeff(k,:),m,m)*Vvis.'; %#ok<NASGU>
    end
end
T.reconstruct=T.reconstruct+toc(phase);
T.total=toc(wall);

stages=[T.setup+T.data,T.fhat,T.dct,T.projection, ...
        T.assembly,T.linearSolve,T.accum,T.reconstruct];
overhead=T.total-sum(stages);

% Account for loop/dispatcher/timing overhead transparently in setup.
% All stages still originate from ONE measured run.
if overhead< -max(1e-9,1e-5*T.total)
    error('Detailed phase times exceed measured total by %g s.',-overhead);
end
stages(1)=stages(1)+overhead;
T.setup=stages(1);
T.stage=stages;
T.stageNames={'setup_and_spatial_data','fhat_or_symbol_evaluation', ...
    'CGL_transform','RHS_vector_assembly','matrix_assembly', ...
    'linear_system_mldivide','CIM_accumulation', ...
    'solution_reconstruction'};

if abs(sum(T.stage)-T.total)>max(1e-9,1e-9*T.total)
    error('The sum of detailed timing components differs from wall time.');
end

out.M=M;
out.N=N;
out.dim=dim;
out.Q=numel(times);
out.dof=ndof;
out.times=times;
out.coeff=coeff;
out.T=T;
% Exact solution errors are computed only ONCE per timing configuration
% by c9_source_error.m, outside all repeated timed trials.
out.errorL2=NaN;

% Coefficient storage is exactly 8*Q*(M-1)^dim bytes (real double).
% If all N complex contour-node solutions were stored, they would require
% an additional 16*N*(M-1)^dim bytes.
out.coefficientBytes=8*out.Q*ndof;
out.allNodeExtraBytes=16*N*ndof;
end

function [nz,dz]=laplace_terms(z,p)
nz=1+p.a*z.^p.alpha;
dz=1+p.b*z.^p.beta;
for k=1:numel(p.aa)
    nz=nz+p.aa(k)*z.^p.ax(k);
end
for k=1:numel(p.bb)
    dz=dz+p.bb(k)*z.^p.bx(k);
end
end

function B=make_basis(M)
% Original CGL-to-Legendre load projection and exactly banded mass matrix.
[x,w]=gauss_legendre(M+1);
P=legendre_dirichlet(x,M);
Cheb=chebvander(x,M);
B.project=P.'*(bsxfun(@times,w,Cheb));
B.nodes=cos(pi*(0:M)'/M);

j=(0:M-2)';
c=1./sqrt(4*j+6);
diagTerms=c.^2.*(2./(2*j+1)+2./(2*j+5));
m=M-1;
main=(1:m)';
mass=sparse(main,main,diagTerms,m,m);
if m>=3
    k=(1:m-2)';
    off=-c(k).*c(k+2).*2./(2*j(k)+5);
    mass=mass+sparse(k,k+2,off,m,m)+sparse(k+2,k,off,m,m);
end
B.mass=mass;
B.stiff=speye(m);
end

function P=legendre_all(x,M)
x=x(:);
P=zeros(numel(x),M+1);
P(:,1)=1;
if M>=1
    P(:,2)=x;
end
for k=1:M-1
    P(:,k+2)=((2*k+1)*x.*P(:,k+1)-k*P(:,k))/(k+1);
end
end

function V=legendre_dirichlet(x,M)
P=legendre_all(x,M);
j=0:M-2;
V=bsxfun(@rdivide,P(:,1:M-1)-P(:,3:M+1),sqrt(4*j+6));
end

function Cheb=chebvander(x,M)
x=x(:);
Cheb=zeros(numel(x),M+1);
Cheb(:,1)=1;
if M>=1
    Cheb(:,2)=x;
end
for k=1:M-1
    Cheb(:,k+2)=2*x.*Cheb(:,k+1)-Cheb(:,k);
end
end

function [x,w]=gauss_legendre(n)
j=(1:n-1)';
off=j./sqrt(4*j.^2-1);
Jac=diag(off,1)+diag(off,-1);
[V,D]=eig(Jac);
[x,ind]=sort(diag(D));
V=V(:,ind);
w=2*(V(1,:).').^2;
end

function C=cosine_dense(M)
% The original direct DCT-I: twice-weighted INTERIOR CGL sample columns,
% with half-weighting of the first/last OUTPUT coefficients.
k=0:M;
C=cos((pi/M)*(k.'*k));
C(:,2:end-1)=2*C(:,2:end-1);
end

function c=cheb_fft1(y)
M=size(y,1)-1;
extension=[y;y(end-1:-1:2,:)];
values=fft(extension,[],1);
c=values(1:M+1,:)/M;
c([1 end],:)=c([1 end],:)/2;
end

function c=cheb_coeff(y,method,C)
M=size(y,1)-1;
if isvector(y)
    if strcmpi(method,'fft')
        c=cheb_fft1(y(:));
    else
        c=(C*y(:))/M;
        c([1 end])=c([1 end])/2;
    end
else
    if strcmpi(method,'fft')
        c=cheb_fft1(y);
        c=cheb_fft1(c.').';
    else
        c=(C*y*C.')/(M*M);
        c([1 end],:)=c([1 end],:)/2;
        c(:,[1 end])=c(:,[1 end])/2;
    end
end
end

function [b,dctTime,projectionTime]=modal_load(B,y,method,C)
% Separates the FFT/direct DCT-I from Chebyshev-to-Legendre projection.
clock=tic;
ch=cheb_coeff(y,method,C);
dctTime=toc(clock);
clock=tic;
if isvector(ch)
    b=B.project*ch(:);
else
    b=B.project*ch*B.project.';
    b=b(:);
end
projectionTime=toc(clock);
end

function [v,vpp]=shape1(x,name)
if strcmpi(name,'analytic')
    h=1.1+sin(pi*x);
    dh=pi*cos(pi*x);
    ddh=-pi^2*sin(pi*x);
elseif strcmpi(name,'finite')
    u=x-.173; q=2.6;
    h=.7+abs(u).^q;
    dh=q*sign(u).*abs(u).^(q-1);
    ddh=q*(q-1)*abs(u).^(q-2);
else
    error('Unknown manufactured shape: %s',name);
end
v=(1-x.^2).*h;
vpp=-2*h-4*x.*dh+(1-x.^2).*ddh;
end

function v=six_gaussians(X,Y)
v=zeros(size(X));
for k=0:5
    cx=.8*cos(k*pi/3);
    cy=.8*sin(k*pi/3);
    v=v+exp(-64*((X-cx).^2+(Y-cy).^2));
end
v=(1-X.^2).*(1-Y.^2).*v;
end

function [z,zp,tau]=hyperbola_nodes(N,t0,Lambda)
% Unchanged contour selection from uploaded comment9_verify.m.
nu=.5;
theta=pi/4;
d=pi/10;
delta=pi/8;
assert(d>0 && d<min(theta,pi/2-theta-delta), ...
    'Contour strip is not admissible.');
rho=linspace(.05,.95,181);
av=acosh(Lambda./((1-nu)*(1-rho)*sin(theta-d)));
epsN=exp(-2*pi*d*N./av);
obj=eps*epsN.^(rho-1)+epsN.^rho./(1-epsN);
[~,j]=min(obj);
r=rho(j);a=av(j);
tau=a/N;
mu=2*pi*d*N*(1-r)/(t0*Lambda*a);
phi=((0:N-1)+.5)*tau;
z=mu*(1+sin(1i*phi-theta));
zp=1i*mu*cos(1i*phi-theta);
end

function self_test()
% Checks no unphysical FFT/direct discrepancy is introduced.
M=16;
rng(17);
v=randn(M+1,1)+1i*randn(M+1,1);
C=cosine_dense(M);
cfft=cheb_coeff(v,'fft',[]);
cdirect=cheb_coeff(v,'direct',C);
a=norm(cfft-cdirect);
assert(a<1e-11,'DCT-I verification failed.');
A=randn(M+1,M+1);
fa=cheb_coeff(A,'fft',[]);
da=cheb_coeff(A,'direct',C);
b=norm(fa-da,'fro');
assert(b<1e-10,'2D DCT-I verification failed.');
B=make_basis(M);
boundary=max(max(abs(legendre_dirichlet([-1;1],M))));
assert(boundary<1e-13);
% Validate the positive-half contour against L^{-1}{1/z}=1.
[z,zp,tau]=hyperbola_nodes(40,.25,3);
test=(tau/pi)*imag(sum(exp(.5*z).*zp./z));
assert(abs(test-1)<1e-9,'Contour inversion check failed.');
fprintf('DCT-I FFT/direct 1D discrepancy: %.3e\n',a);
fprintf('DCT-I FFT/direct 2D discrepancy: %.3e\n',b);
fprintf('Dirichlet endpoint check: %.3e\n',boundary);
fprintf('Contour inversion of 1/z: %.12f\n',test);
fprintf('1D DOF M-1=%d; mass nnz=%d\n',M-1,nnz(B.mass));

% Direct check of MATLAB's mldivide algebraic residual.
A=sparse([4+1i,1;1,3+2i]);
rhs=[1+2i;2-1i];
sol=A\rhs;
relativeResidual=norm(A*sol-rhs)/norm(rhs);
assert(relativeResidual<1e-12, ...
    'MATLAB matrix left division failed the residual check.');
fprintf('A\\b linear-solve relative residual: %.3e\n',relativeResidual);

% Compare complete FFT and non-FFT CIM--CLG solves using A\b.
p.alpha=.50;p.beta=.35;p.gamma=.45;
p.a=1;p.b=.5;p.ax=[];p.aa=[];p.bx=[];p.bb=[];
p.kappa=.25;p.t0=.25;p.Lambda=3;
u1=c9_source_solver(10,24,.5,1,'analytic',p,'fft');
u2=c9_source_solver(10,24,.5,1,'analytic',p,'direct');
relative1D=norm(u1.coeff-u2.coeff)/max(norm(u1.coeff),eps);
assert(relative1D<1e-9, ...
    'FFT and non-FFT one-dimensional solves disagree.');
fprintf('1D FFT/direct full-solver difference: %.3e\n',relative1D);

p.ax=[.10 .18];p.aa=[.35 .20];
p.bx=[.05 .10];p.bb=[.25 .15];
u3=c9_source_solver(8,20,.5,2,'sixpeak',p,'fft');
u4=c9_source_solver(8,20,.5,2,'sixpeak',p,'direct');
relative2D=norm(u3.coeff-u4.coeff)/max(norm(u3.coeff),eps);
assert(relative2D<1e-9, ...
    'FFT and non-FFT two-dimensional solves disagree.');
fprintf('2D FFT/direct full-solver difference: %.3e\n',relative2D);
end
