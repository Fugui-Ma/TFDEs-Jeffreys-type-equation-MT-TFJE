function [modes,cim] = Jeffreys2D_CIM_ModalTimes_KJ2(sp,p,N,times)
% Compute ALL observation times on the user's ORIGINAL hyperbolic CIM contour.
% modes(:,it) is the vectorization of the mass-orthonormal modal expansion:
%   V'*Mm*V=I,  u_coef=V*reshape(modes(:,it),m,m)*V'.
% Therefore, for fixed M, the genuine 2D L2 difference between solutions
% with different N is exactly norm(modes(:,i)-modes_ref(:,i),2).
% No interpolation, artificial reference M, or approximate 2D quadrature.
validateattributes(N,{'numeric'},{'scalar','integer','positive'});
times=times(:).';
if isempty(times) || any(~isfinite(times)) || any(times<=0)
    error('times must contain finite positive observation times.');
end
T=p.t0*p.B_gamma;
if min(times)<p.t0*(1-1e-12) || max(times)>T*(1+1e-12)
    error('All observation times must lie in [t0,t0*B_gamma].');
end
alpha_star=pi/4;
d_bar=pi/8;
d=min(alpha_star,pi/2-d_bar-alpha_star);
[theta,a_theta]=H_optimal_param_out(N,alpha_star,d,p.B_gamma);
tau=a_theta/N;
mu=2*pi*d*N*(1-theta)/(T*a_theta);
k=0:N-1;
s=(k+0.5)*tau;
z=mu*(1+sin(1i*s-alpha_star));
z_prime=1i*mu*cos(1i*s-alpha_star);
Lsum=bsxfun(@plus,sp.lambda,sp.lambda.');
m=sp.M-1;
H=complex(zeros(m*m,N));
for j=1:N
    [P,Q]=Jeffreys2D_PQ(z(j),p);
    numerator=z(j)^(p.Gamma-1)*P;
    denominator=z(j)^p.Gamma*P+Q*Lsum;
    R=(numerator./denominator)*z_prime(j);
    H(:,j)=R(:);
end
% Modal coefficients at all times in one matrix product.
response=(tau/pi)*imag(H*exp(z(:)*times));
modes=bsxfun(@times,sp.Fmode(:),response);
cim=struct('theta',theta,'a_theta',a_theta,'tau',tau,'mu',mu, ...
   'alpha_star',alpha_star,'d_bar',d_bar,'d',d,'z',z,'z_prime',z_prime,...
   't0',p.t0,'B_gamma',p.B_gamma,'T',T,'N',N);
end
