function UUU=Code_2D_solver_zero_Initial_Smooth_CIM(a,b,Alpha,Beta,Gamma,M,Nmax,current_time)
%**************************************************************************
%**    Program： By using the CGL-method to solve the 2D Case :
%                with zero Dirichlet bundary condition
%            
%Explainations:
%   this expriment results showing in  Table 3 and 4
%***Author：
%           F.G. MA
%*****Finishi Data：
%           2022-11-10
%**************************************************************************
% clear variables;
% clear all; close all; clc  
% %parameters
% chi=1; 
% Alpha=0.4; Beta=0.4; Gamma=0.4;
% %Input to this phase is the number of values n
for N=1:Nmax; 
N
% N=50; M=12; current_time=0.5; 
%for iterate=1:length(Mkk)
% set up collocation points x_j: x(j)=cos(pi*j/M), 0<=j<=M
%x=cos(pi*(0:M)'/M); % Chebychev-Gauss-Labotta points
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                   Realize the integal contour of CIM
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
t0=0.01; B_gamma=500; % deal with the time level

T=t0*B_gamma;

% parameters of the actural case for CIM
%d_bar=0.7901; alpha_star=0.6767;
d_bar=pi/8; alpha_star=pi/4;
d=min(alpha_star,pi/2-d_bar-alpha_star);

tic % 计时开始

[theta,a_theta]=H_optimal_param_out(N,alpha_star,d,B_gamma);

%the parameter h of hyperbolic cantour  
tau = a_theta./N;          

%the parametere mu of the hyper contour
mu = 2*pi*d*N*(1-theta)./(T*a_theta);  
clear k
% Hyperbolic contour
for k=0:N-1
        z(k+1) = mu*(1+sin(1i*(k+1/2)*tau-alpha_star));
    z_bar(k+1) = mu*cos(1i*(k+1/2)*tau-alpha_star)*1i;
end
clear k
%++++++++++++++++++++++++++++++++++++++++++++++++++
%    Assemble stiffness matrix, mass matrix
%+++++++++++++++++++++++++++++++++++++++++++++++++++
% stiffness matrix
kk = (0:M-2)';
S = spdiags(4*kk+6,0,M-1,M-1);

% mass matrix (pentadiagonal matrix) 
v1 = -2./(2*kk+5);
Mm = spdiags([v1,2./(2*kk+1)+2./(2*kk+5),[0;0;v1(1:M-3)]],...
             [-2,0,2],(M-1),(M-1));
         
% Kronker integral
A = (kron(S,Mm)+kron(Mm,S));
Mk = kron(Mm,Mm);

%figure(1), clf, spy(A), drawnow

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Evaluate the Legendre coefficients of I^L_Mf(x) from {f(x_j)},j=0,1,...N
%    (backward/forward Chebyshev-Legendre transform)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% set up LGL points and weights
% xx=LGL_points(M)
% w1=LGL_weights(M,xx)
LGLnodes=LGL_points_Weights(M);

[x,y] = meshgrid(LGLnodes,LGLnodes);
   xx = x';
   xx = xx(:);
   yy = y';
   yy = yy(:);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Evaluate the Legendre coefficients of I^L_Mf(x) from {f(x_j)},j=0,1,...N
% (backward/forward Legendre transform: see J. Shena and T. Tang, Science 
%  press, Beijing, 2006.)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% values at the LGL points
    f = Initial_values(xx,yy);
    
    %coef of I_M F w.r.t Legendre polynomials
n1 = M+1;
fcof = zeros(n1,n1); 

%weights of Legendre-Gauss-Labato quadure
w = (2./(M*(M+1))).*(legendreP(M,LGLnodes).^(-2));  
%计算LGL数值积分时对应的求和号外面的权
 
%++++++++++++++++++++++++++++++++++++++++++++++
% Backward legendre transform to get g_hat_{k,j}.
%++++++++++++++++++++++++++++++++++++++++++++++
    for k=1:M+1
        lky = legendreP(k-1,LGLnodes);
        if k == M+1        
           gk = 2./M;
        else
           gk = 2./(2*k-1);  
        end
        for j = 1:M+1
            llx = legendreP(j-1,LGLnodes);
            if j == M+1
               gl = 2./M;   
            else
               gl = 2./(2*j-1); 
            end
               fcof(j,k) = sum(f.*kron(lky.*w,llx.*w))./(gl*gk);
        end
    end


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%    constructor F vector
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
j = (0:M)';
lnorm = 2./(2*j+1); % Coeff. obtained by orthogonality relation 
k1 = 1:M-1;
k2 = 3:n1;
n11 = [(M-1)^2 1];
F = reshape(fcof(k1,k1),n11).*kron(lnorm(k1),lnorm(k1))-...
    reshape(fcof(k1,k2),n11).*kron(lnorm(k2),lnorm(k1))-...
    reshape(fcof(k2,k1),n11).*kron(lnorm(k1),lnorm(k2))+...
    reshape(fcof(k2,k2),n11).*kron(lnorm(k2),lnorm(k2));

% Solve linear system
% 得到相应的5对角阵
tic,
for k=1:N
    % p_hat=(z(k).^(Alpha+Gamma)+z(k).^Gamma).*Mk+...
    %       (1+chi^Beta.*z(k).^Beta).*A;
    % IN=(z(k).^(Alpha+Gamma-1)+z(k).^(Gamma-1)).*F;

     p_hat= (z(k).^(Gamma)+a*z(k).^(Alpha+Gamma))*Mk+...
            (1+b*z(k).^Beta)*A;
      IN = (a*z(k).^(Alpha+Gamma-1)+z(k).^(Gamma-1)).*F;
      
    u_bar{k}=p_hat\IN;
end
toc
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   Forward legendre transforms
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[in3,in4] = meshgrid(0:M-2,0:M-2);
in3 = in3';
in3 = in3(:);
in4 = in4';
in4 = in4(:);
 %uh = zeros((M+1)^2,1);
 for j=1:N
     for k = 1:(M+1)^2
         U{j}(k) = sum(u_bar{j}.*(legendreP(in3,x(k))-legendreP(in3+2,x(k))).*...
                      (legendreP(in4,y(k))-legendreP(in4+2,y(k))));
     end
 end
% At this point, the spatial semi-discrete scheme is complete 
u_h=zeros;
for j=1:N
    u_hatt{j} = exp(z(j)*current_time)*z_bar(j).*U{j};
    u_h = u_h + u_hatt{j};   
end
UU = tau/pi*imag(u_h);

UUU{N} = reshape(UU,M+1,M+1);
end

% lglnodes =LGLnodes;
% 
% 
% % kkk=3;
% % 
% % figure(2)
% % surfc(lglnodes{kkk},lglnodes{kkk},UUU{kkk})
% % shading interp    % interpolate colors across lines and faces
% % colorbar
% % xlabel x, ylabel y, zlabel {Numerical solution}
% % 
% % Interpolate to finer grid and plot:
% [xx,yy] = meshgrid(-1:.01:1,-1:.01:1);
% uuuh = interp2(lglnodes,lglnodes,UUU,xx,yy,'spline');
% figure(3), clf, surfc(xx,yy,uuuh), 
% shading interp    % interpolate colors across lines and faces
% colorbar
% xlabel x, ylabel y, zlabel {Numerical solution}
% 
% figure(4)
% contourf(xx,yy,uuuh,20)
% shading interp    % interpolate colors across lines and faces
% c=colorbar;
% c.Label.String = 'Numerical solution $$p(t,x,y)$$';
% xlabel x, ylabel y
% end
% 
% 
% % % Recall the exact solution to comparing with the numerical solutions
% Exact_solu=Exact_solution(current_time,a,xx);
% 
% % L2 error, L^{\infinity}-error
% Err_inf(1,N) = norm(UU'-Exact_solu,inf);
% Err_L2(1,N) =norm(UU'-Exact_solu,2);
% 
% figure(1)
% plot(xx(2:M),Exact_solu(2:M),'.-',xx(2:M),UU(2:M),'o')
% legend('Exact solution ','CIM-CGL')
% end
% time=toc
% figure(2)
% semilogy(4:4:Nmax,Err_inf(4:4:Nmax),'.--',4:4:Nmax,Err_L2(4:4:Nmax),'.-','markersize',25), 
% %loglog(log(6:2:Mmax),Err_inf,'.:',log(6:2:Mmax),Err_L2,'.-','markersize',12), 
% legend('L^{\infty}-norm','L^2-norm')
% t1=title(['$\alpha=$',num2str(Alpha),',~','$\beta=$',num2str(Beta),',~',...
%     '$\gamma=$',num2str(Gamma),',~','$M=$',num2str(M)]);
% t2=xlabel(' $N$');
% t3=ylabel('$\log_{10}(error_N)$');
% set(t1,'interpreter','latex','FontSize',14);
% set(t2,'interpreter','latex','FontSize',18);
% set(t3,'interpreter','latex','FontSize',18);