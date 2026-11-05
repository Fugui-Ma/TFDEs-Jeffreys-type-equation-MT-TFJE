function [a0,history] = L1_History_Coefficients(t,U,q,n)
%**************************************************************************
% Nonuniform L1 formula for the Caputo derivative at t_n.
%
% U contains the known coefficient vectors [u^0,...,u^(n-1)].
% The result is
%   D_C^q u(t_n) = a0*u^n + history.
%
% For q=1, the backward Euler formula is used.
%**************************************************************************
dt = diff(t(1:n+1));

if abs(q-1) < 1e-12
    a0 = 1/dt(n);
    history = -a0*U(:,n);
    return;
end

tn = t(n+1);
Aq = zeros(n,1);

for k = 1:n
    Aq(k) = ((tn-t(k))^(1-q) ...
            -(tn-t(k+1))^(1-q)) ...
            /(gamma(2-q)*(t(k+1)-t(k)));
end

a0 = Aq(n);
history = -a0*U(:,n);

for k = 1:n-1
    history = history + Aq(k)*(U(:,k+1)-U(:,k));
end
end
