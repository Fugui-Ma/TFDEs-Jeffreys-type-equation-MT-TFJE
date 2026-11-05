function f = source_term_time_L1(kappa,a,b,Alpha,Beta,Gamma,t,x)
%**************************************************************************
% Time-domain source term for the 1D manufactured solution
% u(t,x)=t^kappa*(sin(pi*x)+1.1)*(1-x^2).
%
% It is the inverse Laplace transform of the source_term.m expression
% supplied by the user. The fractional derivatives are Caputo derivatives.
%**************************************************************************
v = (sin(pi*x)+1.1).*(1-x.^2);

v_xx = -sin(pi*x).*(pi^2*(1-x.^2)+2) ...
       -2.2-4*pi*x.*cos(pi*x);

f = gamma(kappa+1)*( ...
      v.*( t^(kappa-Gamma)/gamma(kappa+1-Gamma) ...
         + a*t^(kappa-Alpha-Gamma)/gamma(kappa+1-Alpha-Gamma)) ...
    - v_xx.*(t^kappa/gamma(kappa+1) ...
         + b*t^(kappa-Beta)/gamma(kappa+1-Beta)) );
end
