function [L2_error,time] = CPU_time_cost_CIM_CLG(...
    kappa,a,b,Alpha,Beta,Gamma,N,M,current_time,T,B_gamma)
%**************************************************************************
% CIM-CLG solver for the 1D manufactured-solution problem.
% The spatial CLG basis is phi_k(x)=L_k(x)-L_{k+2}(x).
% The source term and contour parameters use the uploaded original routines.
%**************************************************************************
time_start = tic;

%--------------------------------------------------------------
% Realize the hyperbolic contour (same formula as the original code).
%--------------------------------------------------------------
d_bar = pi/8;
alpha_star = pi/4;
d = min(alpha_star,pi/2-d_bar-alpha_star);

[theta,a_theta] = H_optimal_param_out(N,alpha_star,d,B_gamma);

tau = a_theta/N;
mu = 2*pi*d*N*(1-theta)/(T*a_theta);

z = zeros(1,N);
z_bar = zeros(1,N);
for k = 0:N-1
    z(k+1) = mu*(1+sin(1i*(k+1/2)*tau-alpha_star));
    z_bar(k+1) = mu*cos(1i*(k+1/2)*tau-alpha_star)*1i;
end

%--------------------------------------------------------------
% Spatial CLG matrices and LGL points.
%--------------------------------------------------------------
xx = LGL_points_Weights(M);
[Mm,S] = CLG_Matrices(M);

%--------------------------------------------------------------
% Laplace transform of the source term: original source_term.m.
%--------------------------------------------------------------
f = source_term(kappa,a,b,Alpha,Beta,Gamma,xx,z);

F = zeros(M-1,N);
for i = 1:N
    F(:,i) = CLG_Load_Vector(f(:,i),xx,M);
end

%--------------------------------------------------------------
% Solve the M-1 dimensional complex linear systems.
% The matrices are complex symmetric but generally non-Hermitian.
% Use an explicit general LU decomposition.
%--------------------------------------------------------------
u_bar = complex(zeros(M-1,N));

for k = 1:N
    p_hat = (z(k)^Gamma+a*z(k)^(Alpha+Gamma))*Mm ...
          + (1+b*z(k)^Beta)*S;

    D = decomposition(p_hat,'lu');
    u_bar(:,k) = D\F(:,k);
end

%--------------------------------------------------------------
% CIM accumulation in the CLG coefficient vector.
%--------------------------------------------------------------
u_coef = complex(zeros(M-1,1));
for j = 1:N
    u_coef = u_coef + exp(z(j)*current_time)*z_bar(j)*u_bar(:,j);
end
u_coef = tau/pi*imag(u_coef);

%--------------------------------------------------------------
% L2 error, evaluated independently of the solution nodes.
%--------------------------------------------------------------
L2_error = CLG_L2_Error(u_coef,M,kappa,current_time);

% CPU time includes the error evaluation, as in the original program.
time = toc(time_start);
end
