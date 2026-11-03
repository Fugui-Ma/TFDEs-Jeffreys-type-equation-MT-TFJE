function [Err_inf,Err_L2]=Space_error_order_CGL(kappa,a,b,Alpha,Beta,Gamma,N,M,current_time)
%**************************************************************************
%**    Program： By using the CGL-method to solve the 1D Case :
%                with zero Dirichlet bundary condition
%            
%Explainations:
%   this expriment results showing in  Table 3 and 4
%***Author：
%           F.G. MA
%*****Finishi Data：
%           2022-8-13
%**************************************************************************
%clear all; close all; clc  

% parameters
% a=3/2; chi=100; 
% Alpha=0.5; Beta=0.5; Gamma=0.5;

% Input to this phase is the number of values n
%N=100; 
current_time=current_time; 
% Mmax=26; 
% for M=6:2:Mmax
M
% set up collocation points x_j: x(j)=cos(pi*j/M), 0<=j<=M
%x=cos(pi*(0:M)'/M); % Chebychev-Gauss-Labotta points
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                   Realize the integal contour of CIM
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
t0=0.1; B_gamma=10; % deal with the time level

T=t0*B_gamma;

% parameters of the actural case for CIM
%d_bar=0.7901; alpha_star=0.6767;
d_bar=pi/8; alpha_star=pi/4;
d=min(alpha_star,pi/2-d_bar-alpha_star);

tic % 计时开始

[theta,a_theta]=H_optimal_param_out(N,alpha_star,d,B_gamma);

%the parameter h of hyperbolic cantour  
tau = a_theta/N;          

%the parametere mu of the hyper contour
mu = 2*pi*d*N*(1-theta)/(T*a_theta);  

% Hyperbolic contour
for k=0:N-1
        z(k+1) = mu*(1+sin(1i*(k+1/2)*tau-alpha_star));
    z_bar(k+1) = mu*cos(1i*(k+1/2)*tau-alpha_star)*1i;
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Pre-computation: Compute LGL points, a_k,b_k,S_kk
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%[x]=LGL_points_Weights(M);  %x, LGL points. w: weights. 
xx=LGL_points_Weights(M);

% parameters for Legendre basis
ak=zeros(1,M-1); 
bk=-1.*ones(1,M-1);

%++++++++++++++++++++++++++++++++++++++++++++++++++
%    Assemble stiffness matrix, mass matrix
%+++++++++++++++++++++++++++++++++++++++++++++++++++
% stiffness matrix
S=zeros(M-1,M-1);
for k=0:M-2
    for j=0:M-2
        if k==j
           S(k+1,j+1)=-(4*k+6)*bk(k+1);
        end
    end
end

% mass matrix (pentadiagonalmatrix) 
Mm=zeros(M-1,M-1);
for j=0:M-2
    for k=0:M-2
         if j==k
             Mm(j+1,k+1) = 2./(2*k+1)+ak(k+1)^2*2./(2*k+3)+...
                           bk(k+1)^2*2./(2*k+5);
         elseif j==k+1
             Mm(j+1,k+1) = ak(k+1)*2./(2*k+3)+ak(k+2)*2./(2*k+5);
             Mm(k+1,j+1) = Mm(j+1,k+1);
         elseif j==k+2
             Mm(j+1,k+1) = bk(k+1)*2./(2*k+5);
             Mm(k+1,j+1) = Mm(j+1,k+1);
         else
             Mm(j+1,k+1) = 0;
         end
    end
end
%figure, clf, spy(Mm),drawnow   % Show the structure of the Mm.
% Compute the eigenvalues of Mm
% [E_h,lambda]=eig(Mm);
% % Construct the Fourier-like basis
%  for k=1:M-1
%      for j=0:
%      phi(k)
%  end

% Mass matrix generated with the Fourier-like basis
% MM=zeros(M-1,M-1);
% for k=0:M-2
%     for j=0:M-2
%         if k==j
%            MM(k+1,j+1)=lambda(k+1);
%         end
%     end
% end 
%MM=lambda

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Evaluate the Legendre coefficients of I^L_Mf(x) from {f(x_j)},j=0,1,...N
% (backward/forward Legendre transform: see J. Shena and T. Tang, Science 
%  press, Beijing, 2006.)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% set up the coefficients c_k: c(0)=2,c(M)=2, c(j)=1,1<=j<=N-1
%c=[2;ones(M-1,1);2]';
% values at the LGL points
[f]=source_term(kappa,a,b,Alpha,Beta,Gamma,xx,z);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Backward Legendre transform to get g_hat_k. Cf. P.112
% Recall the Legendre polynomial of degree n 
%   y=Legendre_polynomial(n,x)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
g_hat=zeros(M+1,N);
for i=1:N
    for k=0:M-1
        temp=zeros;
        for j=0:M
           % Legendre polynomial of degree k
            L_k=Legendre_polynomial(k,xx(j+1));
           % Legendre polynomial of degree  M
            L_M=Legendre_polynomial(M,xx(j+1));
           % Backward Legendre transform
            temp=temp+f(j+1,i).*L_k./(L_M).^2;
        end
        g_hat(k+1,i)=(2*k+1)./(M*(M+1)).*temp;
    end
end
clear temp
%g_hat
for i=1:N
    temp=zeros;
        for j=0:M
           % Legendre polynomial of degree  M
            L_M=Legendre_polynomial(M,xx(j+1));
           % Backward Legendre transform
            temp=temp+f(j+1,i).*1.0./L_M;
        end
        g_hat(M+1,i)=1.0./(M+1).*temp;
end
clear temp
%g_hat
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Evaluate f_bar from f_k=(\Sum_{j=0}^M*g_hat_j*L_j(x),\phi_j(x))
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear F
for k=0:M-2
    F(k+1,:)= g_hat(k+1,:)/(k+1/2)+ak(k+1)*g_hat(k+2,:)/(k+3/2)+...
        bk(k+1)*g_hat(k+3,:)/(k+5/2);
end
% Solve linear system

for k=0:N-1
    %p_hat=(z(k+1).^(Alpha+Gamma)+z(k+1).^(Gamma))*Mm+...
    %           (1+chi^Beta*z(k+1).^Beta)*S;
    p_hat=(z(k+1).^(Gamma)+a*z(k+1).^(Alpha+Gamma))*Mm+...
            (1+b*z(k+1).^Beta)*S;
    u_bar{k+1}=p_hat\F(:,k+1);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Evaluate g_k from \Sum_{j=0}^{M-2}*u_bar_j*\phi_j(x) 
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for j=1:N
     g{j}(1)=u_bar{j}(1); g{j}(2)=u_bar{j}(2)+ak(1)*u_bar{j}(1);
        for  k=2:M-2
         g{j}(k+1)=u_bar{j}(k+1)+ak(k)*u_bar{j}(k)+bk(k-1)*u_bar{j}(k-1);
        end
     g{j}(M)=ak(M-1)*u_bar{j}(M-1)+bk(M-2)*u_bar{j}(M-2);   
     g{j}(M+1)=bk(M-1)*u_bar{j}(M-1);
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Forward legendre transform
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for i=1:N
    for k=0:M
    temp=zeros;
        for j=0:M
            L_k=Legendre_polynomial(j,xx(k+1));
            temp=temp+g{i}(j+1)*L_k;
        end
    U{i}(k+1)=temp;
    end
end
% At this point, the spatial semi-discrete scheme is complete 
u_h=zeros;
for j=1:N
    u_hatt{j} = exp(z(j)*current_time)*z_bar(j).*U{j};
      u_h = u_h + u_hatt{j};   
end
UU = tau/pi*imag(u_h);

% Recall the exact solution to comparing with the numerical solutions
Exact_solu=Exact_solution(current_time,kappa,xx);

% L2 error, L^{\infinity}-error
% Evaluate errors independently of the M+1 solution nodes.
M_error=max(128,2*M+2);
[x_error,w_error]=LGL_points_Weights(M_error);
g_time=zeros(M+1,1);
for j=1:N
    g_time=g_time+exp(z(j)*current_time)*z_bar(j).*g{j}(:);
end
g_time=tau/pi*imag(g_time);
UU_error=zeros(size(x_error));
for k=0:M
    UU_error=UU_error+g_time(k+1)*Legendre_polynomial(k,x_error);
end
Exact_solu_error=Exact_solution(current_time,kappa,x_error);
L2_error=sqrt(sum(w_error.*abs(UU_error-Exact_solu_error).^2));
% Dense-grid approximation of the continuous maximum norm.
x_inf=linspace(-1,1,20001)';
UU_inf=zeros(size(x_inf));
for k=0:M
    UU_inf=UU_inf+g_time(k+1)*Legendre_polynomial(k,x_inf);
end
Exact_solu_inf=Exact_solution(current_time,kappa,x_inf);
Err_inf = max(abs(UU_inf-Exact_solu_inf));
Err_L2 = L2_error;

figure(1)
plot(xx(2:M),Exact_solu(2:M),'.-',xx(2:M),UU(2:M),'o')
legend('Exact solution ','CIM-CGL')

return
end
%end
% time=toc
% figure(2)
% semilogy(6:2:Mmax,Err_inf,'.--',6:2:Mmax,Err_L2,'.-','markersize',25), 
% %loglog(log(6:2:Mmax),Err_inf,'.:',log(6:2:Mmax),Err_L2,'.-','markersize',12), 
% legend('L^{\infty}-norm','L^2-norm')
% t1=title(['$\alpha=$',num2str(Alpha),',~','$\beta=$',num2str(Beta),',~',...
%     '$\gamma=$',num2str(Gamma),',~','$N=$',num2str(N)]);
% t2=xlabel('polynomial degree $M$');
% t3=ylabel('$\log_{10}(error_M)$');
% set(t1,'interpreter','latex','FontSize',14);
% set(t2,'interpreter','latex','FontSize',18);
% set(t3,'interpreter','latex','FontSize',18);


