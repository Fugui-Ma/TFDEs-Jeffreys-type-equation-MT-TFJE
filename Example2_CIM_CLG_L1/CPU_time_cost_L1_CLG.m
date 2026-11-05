function [L2_error,time] = CPU_time_cost_L1_CLG(...
    kappa,a,b,Alpha,Beta,Gamma,Ntilde,M,current_time,grade)
%**************************************************************************
% L1-CLG solver for the 1D manufactured-solution problem.
% The spatial CLG basis is exactly the same as in CIM_time_cost_CIM_CLG.
% A graded time mesh is used to handle the initial singularity t^kappa.
%**************************************************************************
time_start = tic;

%--------------------------------------------------------------
% Spatial CLG matrices and LGL points.
%--------------------------------------------------------------
xx = LGL_points_Weights(M);
[Mm,S] = CLG_Matrices(M);

% Graded temporal mesh on [0,current_time].
t = current_time*((0:Ntilde)/Ntilde).^grade;

% U(:,n+1) stores the CLG coefficients at t_n.
% The initial value is zero, consistent with Exact_solution.m.
U = zeros(M-1,Ntilde+1);

%--------------------------------------------------------------
% L1 time stepping.
%--------------------------------------------------------------
for n = 1:Ntilde

    [cGamma,hGamma] = L1_History_Coefficients(t,U(:,1:n),Gamma,n);

    [cAlphaGamma,hAlphaGamma] = ...
        L1_History_Coefficients(t,U(:,1:n),Alpha+Gamma,n);

    [cBeta,hBeta] = L1_History_Coefficients(t,U(:,1:n),Beta,n);

    % Source term in the time domain, obtained from source_term.m
    % by inverse Laplace transformation for the manufactured solution.
    f_n = source_term_time_L1(kappa,a,b,Alpha,Beta,Gamma,t(n+1),xx);
    F_n = CLG_Load_Vector(f_n,xx,M);

    % Fully discrete CLG system at t_n.
    A = (cGamma+a*cAlphaGamma)*Mm + (1+b*cBeta)*S;

    rhs = F_n ...
        - Mm*(hGamma+a*hAlphaGamma) ...
        - b*S*hBeta;

    % This matrix is real symmetric positive definite.
    U(:,n+1) = A\rhs;
end

%--------------------------------------------------------------
% L2 error at the final time.
%--------------------------------------------------------------
u_coef = U(:,end);
L2_error = CLG_L2_Error(u_coef,M,kappa,current_time);

% CPU time includes the error evaluation.
time = toc(time_start);
end
