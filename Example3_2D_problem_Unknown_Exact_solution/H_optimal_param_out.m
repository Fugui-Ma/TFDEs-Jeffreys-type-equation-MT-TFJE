function [theta,a_theta] = H_optimal_param_out(n,alpha,d,B_gamma)
% Retained from the user's original MATLAB helper (numeric operation only).
epsilon=2.22*10^(-16);
theta=linspace(0,1,202);
theta=theta(2:201);
a_theta=acosh(B_gamma./((1-theta)*sin(alpha-d)));
epsilon_theta=exp(-2*pi*d./a_theta.*n);
const=theta.*epsilon+(1-theta).*epsilon_theta./(1-epsilon_theta);
[~,indx]=min(const);
theta=theta(indx);
a_theta=a_theta(indx);
end
