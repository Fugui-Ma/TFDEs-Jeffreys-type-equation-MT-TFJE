function result = Solve_Jeffreys2D_KJ2(sp,cfg)
% Equivalent to solving at each CIM point:
% [z^gamma P(z) kron(Mm,Mm) + Q(z)(kron(S,Mm)+kron(Mm,S))] uhat
%   = z^(gamma-1) P(z) F(:).
% The tensor generalized eigensolver is an exact algebraic reduction.
validateattributes(cfg.N,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.t,{'numeric'},{'scalar','positive'});
assert(cfg.K==2 && cfg.J==2,'This experiment is configured for K=J=2.');
assert(numel(cfg.alpha_terms)==cfg.K && numel(cfg.a_terms)==cfg.K);
assert(numel(cfg.beta_terms)==cfg.J && numel(cfg.b_terms)==cfg.J);
assert(all(cfg.alpha_terms>0 & cfg.alpha_terms<cfg.Alpha), ...
       'Require 0 < alpha_k < alpha.');
assert(all(cfg.beta_terms>0 & cfg.beta_terms<cfg.Beta), ...
       'Require 0 < beta_j < beta.');
assert(cfg.Alpha>0 && cfg.Alpha<=1 && cfg.Beta>0 && cfg.Beta<=1 ...
       && cfg.Gamma>0 && cfg.Gamma<=1);
assert(cfg.a>0 && cfg.b>0 && all(cfg.a_terms>=0) && all(cfg.b_terms>=0));

T = cfg.t0*cfg.B_gamma;
if cfg.t<cfg.t0 || cfg.t>T
    warning('Observation time t=%g lies outside recommended CIM window [%g,%g].', ...
             cfg.t,cfg.t0,T);
end
alpha_star=pi/4;
d_bar=pi/8;
d = min(alpha_star,pi/2-d_bar-alpha_star);
[theta,a_theta] = H_optimal_param_out(cfg.N,alpha_star,d,cfg.B_gamma);
tau = a_theta/cfg.N;
mu = 2*pi*d*cfg.N*(1-theta)/(T*a_theta);
k = 0:cfg.N-1;
s = (k+1/2)*tau;
z = mu*(1+sin(1i*s-alpha_star));
z_prime = mu*cos(1i*s-alpha_star)*1i;

lam_sum = bsxfun(@plus,sp.lambda,sp.lambda.');
G = complex(zeros(size(lam_sum)));
for j=1:cfg.N
    [P,Q] = Jeffreys2D_PQ(z(j),cfg);
    denom = z(j)^cfg.Gamma*P + Q*lam_sum;
    numerator = z(j)^(cfg.Gamma-1)*P;
    % Exactly the same upper-branch CIM quadrature tau/pi*Im(sum(...)).
    G = G + exp(z(j)*cfg.t)*z_prime(j)*(numerator./denom);
end
G = (tau/pi)*imag(G);
Ucoef = sp.V * (sp.Fmode.*G) * sp.V.';
U_xy = sp.basis_nodes * Ucoef * sp.basis_nodes.';

% Internally the coefficient matrix rows index x and columns index y.
% For contourf(meshgrid(nodes,nodes),U), rows MUST index y and columns x.
result = struct();
result.U = real(U_xy.');
% Save the exact Galerkin expansion for accurate common-grid comparisons.
% Ucoef(i,j) multiplies phi_i(x)*phi_j(y), with i,j=0,...,M-2.
result.Ucoef = real(Ucoef);
result.nodes = sp.nodes;
result.p0_nodes = sp.p0_nodes;
result.F = sp.F;
result.weights = sp.weights;
result.params = cfg;
result.params.M = sp.M;
result.cim = struct('alpha_star',alpha_star,'d_bar',d_bar,'d',d, ...
    'theta',theta,'a_theta',a_theta,'tau',tau,'mu',mu, ...
    'T',T,'z',z,'z_prime',z_prime);
result.initial_LGL_mass = sp.initial_LGL_mass;
end
