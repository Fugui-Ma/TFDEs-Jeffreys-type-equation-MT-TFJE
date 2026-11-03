function [theta,a_theta]=H_optimal_param_out(n,alpha,d,B_gamma)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Program: compute the optimal parameters (theta) and the convergence rate (const);
%          note that 0< theta <1; we have modified the range of alpha; 
%
% Author:
%         Ma, F.G
%
% Data:
%         2021-6-6
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%clear;clc

%B_gamma=5;

%delta=0;

%alpha = 0.7;

%n=50;

%d = 0.6;
    
epsilon=2.22*10^(-16);        % machine precision

theta=linspace(0,1,202);

theta=theta(2:201);           % discard the end-point;
 
a_theta = acosh(B_gamma./((1-theta)*sin(alpha-d)));

epsilon_theta=exp(-2*pi*d./a_theta*n);

const=epsilon*epsilon_theta.^(theta-1)+epsilon.^theta;

[const,indx]=min(const);

theta=theta(indx);

a_theta=a_theta(indx);



