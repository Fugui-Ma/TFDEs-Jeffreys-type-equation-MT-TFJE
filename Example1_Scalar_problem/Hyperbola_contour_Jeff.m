function [Err,mu]=Hyperbola_contour_Jeff(N,t,T,alpha,beta,gama,a,b,d,alpha1,B,u_0)

%Program: Using Hyperbolic Contour integral method solve TF-TFJE
%      
%           
%Source£º 
%
% Author£º Ma.F.G
%   Data£º
%          2019-10-14
%  
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%clear all; clc
 
%  N = 150;
%  t0 = 1;          % The time
  B_gamma = B; lambda=1.5;
%  t = 1;    %The time interval of the integral.   
  
  
%  alpha = 0.7
  
%  d=0.6;
%  alpha1=0.7;  %angle in the 
%  u_0 = 1;   % initial value
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% caculus the optimal parameters of the Hyperbolic contour  
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   delta = 0;       % used for define the slope of an asymptotic line
%   theta=0.9;   %theta include in the interval (0,1) 

for j = 1:N
    j
    
    u_h = 0;
    
    [theta,a_theta]=H_optimal_param_out(j,alpha1,d,B_gamma);   % the 1 represents the time interval
    
    h1 = a_theta/j;                                     %the parameter h of hyperbolic cantour
    
    mu = 2*pi*d*j*(1-theta)/(T*a_theta);  %the parametere mu of the hyper contour


     for k = 0:j-1
         
%*************************************************************************************         
% The Hyperbolic cantour integral method solving the Fynman-Kac in two internal tates
%**************************************************************************

  z = mu*(1+sin(1i*(k+1/2)*h1-alpha1));
  z_bar = mu*cos(1i*(k+1/2)*h1-alpha1)*1i;
   
  
  %f_hat = (z+z^alpha/Gamma(1-alpha)+1)*(z^(-1)+3*sqrt(pi)*z^(-3))-(1+z^(alpha-1)/Gamma(1-alpha))*u_0;
  %u_hat_phi = (z+z^alpha/Gamma(1-alpha)+1)^(-1)*(f_hat+(1+z^(alpha-1)/Gamma(1-alpha))*u_0);
     
  % f_hat =4*sqrt(pi)*z.^(alpha+gama-4)+4*sqrt(pi)*z.^(gama-4)+...
  %        lambda*(z.^(-1)+4*sqrt(pi)*z.^(-4)+chi^beta*(z.^(beta-1)+...
  %        2/3*gamma(4)*sqrt(pi)*z.^(beta-4)));
  % 
  % u_hat_phi = (z.^(alpha+gama)+z.^gama+lambda*(1+chi^beta.*z.^beta)).^(-1)*...
  %             (f_hat+(z.^(alpha+gama-1)+z.^(gama-1))*u_0);

f_hat = lambda*z.^(-1) + lambda*gamma(9/5)*z.^(-9/5) + lambda*b*z.^(beta-1)+...
        lambda*b*gamma(9/5)*z.^(beta-9/5) + 4/5*gamma(4/5)*z.^(gama-9/5)+....
        a*4/5*gamma(4/5)*z.^(alpha+gama-9/5);

u_hat_phi = (z.^gama*(1+a*z.^(alpha))+lambda*(1+b*z.^(beta))).^(-1)*...
            (f_hat+z.^(gama-1)*(1+a*z.^(alpha))*u_0);
  
u_hat =  exp(z*t)*z_bar*u_hat_phi;
   
% The sum  
   u_h=u_h+u_hat;

    end 
%the solutions  
    u_H = h1/pi*imag(u_h);  
  
    D1(j) = u_H;
    
    %u_t(j)=1+2*sqrt(pi)/3*t^3;
    u_t(j)=1+t.^(4/5);
 
    Err(j)=abs(u_t(j)-D1(j));
  
end

%exact solution

% figure
% semilogy(1:N,Err,'r.-')
% grid on
% legend('Time error for u(t)')
% xlabel('M');ylabel('log_{10}(max(|u_{h}-u(t)|))');%time_HC = toc
% hold on
% % figure(4)
% % plot(1:N,G1_H,'b-*',1:N,G2_H,'r-p');
% % legend('G_{1,t}','G_{2,t}')
%  
%  %axis([1 N -0.2 0.2])
