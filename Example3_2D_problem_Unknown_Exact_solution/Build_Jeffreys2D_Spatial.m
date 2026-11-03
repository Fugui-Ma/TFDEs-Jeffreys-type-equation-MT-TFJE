function sp = Build_Jeffreys2D_Spatial(M,initial_function)
% SAME Shen Legendre-Galerkin matrices and Legendre coefficient load as the
% supplied Code_2D_solver.m, rearranged as tensor matrices for efficiency.
% phi_k(x)=P_k(x)-P_{k+2}(x), k=0,...,M-2.
if nargin<2 || isempty(initial_function)
    initial_function = @Initial_values_sixpeak_KJ2;
end
[nodes,w] = LGL_points_Weights(M);
[X,Y] = meshgrid(nodes,nodes);
p0_nodes = initial_function(X,Y);  % rows are y, columns are x
if ~isequal(size(p0_nodes),size(X)) || any(~isfinite(p0_nodes(:)))
    error('Initial function must return finite values on the (M+1)x(M+1) LGL grid.');
end

m = M-1;
L = Legendre_Vandermonde_Jeffreys(nodes,M);
% Original fcof(p,q) = sum_i sum_j p0(x_i,y_j) w_i w_j
%                                  P_p(x_i)P_q(y_j) / (g_p g_q).
% LGL is in increasing order; p0_nodes.' indexes (x,y) as in the original.
Lw = bsxfun(@times,L,w);
moments = Lw.' * p0_nodes.' * Lw;
g = 2./(2*(0:M)'+1);
g(end) = 2/M;       % discrete LGL norm for P_M (as in original code)
fcof = moments ./ (g*g.');

% F uses continuous Legendre norms, also exactly as in the old code.
lnorm = 2./(2*(0:M)'+1);
i = 1:m;
ip2 = i+2;
F = fcof(i,i)     .* (lnorm(i)*lnorm(i).') ...
  - fcof(i,ip2)   .* (lnorm(i)*lnorm(ip2).') ...
  - fcof(ip2,i)   .* (lnorm(ip2)*lnorm(i).') ...
  + fcof(ip2,ip2) .* (lnorm(ip2)*lnorm(ip2).');

% Original 1D stiffness and mass matrices.
kk = (0:M-2)';
S = spdiags(4*kk+6,0,m,m);
Mm = spdiags(2./(2*kk+1)+2./(2*kk+5),0,m,m);
for k = 1:(m-2)
    v = -2/(2*(k-1)+5);
    Mm(k,k+2)=v;
    Mm(k+2,k)=v;
end

% S v_i = lambda_i Mm v_i. Normalize so V'*Mm*V=I.
% This diagonalizes the ORIGINAL Kronecker system, not a new discretization.
[V,D] = eig(full(S),full(Mm));
lambda = real(diag(D));
[lambda,idx] = sort(lambda,'ascend');
V = real(V(:,idx));
for k=1:m
    V(:,k) = V(:,k) / sqrt(V(:,k).'*Mm*V(:,k));
end
check = norm(V.'*Mm*V-eye(m),'fro');
if check > 1e-8
    warning('Mass orthogonality residual: %.3e',check);
end

basis_nodes = L(:,1:m)-L(:,3:m+2);
sp.M = M;
sp.nodes = nodes;
sp.weights = w;
sp.p0_nodes = p0_nodes;
sp.F = F;
sp.S = S;
sp.Mm = Mm;
sp.V = V;
sp.lambda = lambda;
sp.Fmode = V.'*F*V;
sp.basis_nodes = basis_nodes;
sp.initial_function = func2str(initial_function);
sp.initial_LGL_mass = sum(sum(p0_nodes .* (w*w.')));
end
