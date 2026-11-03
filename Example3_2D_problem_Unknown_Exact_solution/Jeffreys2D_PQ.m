function [P,Q] = Jeffreys2D_PQ(z,p)
% Multi-term time-fractional Jeffreys factors, K=J=2.
% Complex powers use MATLAB's principal branch, as in original CIM code.
P = 1 + p.a*z.^p.Alpha;
Q = 1 + p.b*z.^p.Beta;
for k = 1:p.K
    P = P + p.a_terms(k)*z.^p.alpha_terms(k);
end
for j = 1:p.J
    Q = Q + p.b_terms(j)*z.^p.beta_terms(j);
end
end
