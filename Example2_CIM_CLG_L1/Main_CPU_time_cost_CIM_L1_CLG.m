%**************************************************************************
% Main program: CPU-time and temporal-resolution comparison
%              between CIM-CLG and L1-CLG for the 1D problem.
%
% The spatial discretization is the same CLG method in both algorithms.
% Figure 1: M versus N (CIM) and N_tilde (L1).
% Figure 2: M versus CPU time.
%
%***Author：
%           F.G. MA
%*****Finishi Data：
%           2026-11-1
%**************************************************************************
clear all; close all; clc

% Parameters: same default set as Main_CPU_time_cost_FGL.m
kappa = 1/4;
a = 1;
b = 1000;

Alpha = 0.75;
Beta  = 0.15;
Gamma = 0.15;

% Other parameter sets in the original programs:
% Alpha = 0.25; Beta = 0.55; Gamma = 0.75;
% Alpha = 0.50; Beta = 0.35; Gamma = 0.45;

current_time = 0.50;

% CIM contour parameters: T=t0*B_gamma, as in the original code.
t0 = 0.1;
B_gamma = 10;
T = t0*B_gamma;

% Prescribed L2 accuracy
tol = 2.2025*10^(-5);

% Spatial polynomial degrees (same range as the supplied comparison plot)
M = 20:2:30;

% Search limits for the temporal discretization parameters
Nmax = 100;                  % maximum CIM contour nodes
Ntilde_max = 4096;           % maximum L1 time steps
Ntilde_test = 16:16:Ntilde_max;  % coarse bracket; refine to an integer below

% L1 graded mesh: t_n=current_time*(n/Ntilde)^r.
% For u(t)~t^kappa, choose r >= (2-qmax)/kappa.
qmax = max([Gamma, Alpha+Gamma, Beta]);
if Gamma <= 0 || Gamma > 1 || Beta <= 0 || Beta > 1 || ...
        Alpha+Gamma <= 0 || Alpha+Gamma > 1
    error('This L1 code requires 0<Gamma, Beta, Alpha+Gamma<=1.');
end
grade = max(1,ceil((2-qmax)/kappa));

% Repeat the selected computations and report the median CPU time.
nRepeat = 3;

fprintf('\n');
fprintf('CIM-CLG versus L1-CLG: one-dimensional manufactured solution\n');
fprintf('L2 tolerance = %.7e\n',tol);
fprintf('current_time = %.6f, T = %.6f\n',current_time,T);
fprintf('Alpha = %.4f, Beta = %.4f, Gamma = %.4f\n',Alpha,Beta,Gamma);
fprintf('kappa = %.4f, a = %.4f, b = %.4f\n',kappa,a,b);
fprintf('L1 graded mesh exponent r = %d\n\n',grade);

% Warm-up calls. These values are not used in the comparison results.
[~,~] = CPU_time_cost_CIM_CLG(kappa,a,b,Alpha,Beta,Gamma,15,12,...
                              current_time,T,B_gamma);
[~,~] = CPU_time_cost_L1_CLG(kappa,a,b,Alpha,Beta,Gamma,64,12,...
                              current_time,grade);

N_CIM = nan(size(M));
Ntilde_L1 = nan(size(M));
Err_CIM = nan(size(M));
Err_L1 = nan(size(M));
CPU_CIM = nan(size(M));
CPU_L1 = nan(size(M));

for kk = 1:length(M)
    M_now = M(kk);
    fprintf('------------------------------------------------------------\n');
    fprintf('M = %d\n',M_now);

    %--------------------------------------------------------------
    % 1. Find the smallest CIM contour-node number N satisfying Tol.
    %--------------------------------------------------------------
    found_CIM = false;
    for N = 1:Nmax
        [error_now,~] = CPU_time_cost_CIM_CLG(kappa,a,b,Alpha,Beta,...
                         Gamma,N,M_now,current_time,T,B_gamma);
        if error_now <= tol
            N_CIM(kk) = N;
            Err_CIM(kk) = error_now;
            found_CIM = true;
            break;
        end
    end

    if ~found_CIM
        warning('CIM-CLG did not reach Tol for M=%d and N<=%d.',...
                 M_now,Nmax);
    end

    %--------------------------------------------------------------
    % 2. Find the smallest L1 time-step number Ntilde satisfying Tol.
    %    First bracket the threshold using 16-step increments, then test every integer in the
    %    first passing bracket.
    %--------------------------------------------------------------
    found_L1 = false;
    Ntilde_previous = 0;

    for Ntilde_candidate = Ntilde_test
        [error_now,~] = CPU_time_cost_L1_CLG(kappa,a,b,Alpha,Beta,...
                    Gamma,Ntilde_candidate,M_now,current_time,grade);

        if error_now <= tol
            for Ntilde = Ntilde_previous+1:Ntilde_candidate
                [error_refine,~] = CPU_time_cost_L1_CLG(kappa,a,b,...
                    Alpha,Beta,Gamma,Ntilde,M_now,current_time,grade);
                if error_refine <= tol
                    Ntilde_L1(kk) = Ntilde;
                    Err_L1(kk) = error_refine;
                    found_L1 = true;
                    break;
                end
            end
            break;
        end
        Ntilde_previous = Ntilde_candidate;
    end

    if ~found_L1
        warning('L1-CLG did not reach Tol for M=%d and Ntilde<=%d.',...
                 M_now,Ntilde_max);
    end

    %--------------------------------------------------------------
    % 3. Measure CPU time at the selected N and Ntilde.
    %    The subprogram timing includes the L2-error evaluation,
    %    consistent with the timing style in the original program.
    %--------------------------------------------------------------
    if found_CIM
        time_samples = zeros(nRepeat,1);
        for ir = 1:nRepeat
            [Err_CIM(kk),time_samples(ir)] = CPU_time_cost_CIM_CLG(...
                kappa,a,b,Alpha,Beta,Gamma,N_CIM(kk),M_now,...
                current_time,T,B_gamma);
        end
        CPU_CIM(kk) = median(time_samples);
    end

    if found_L1
        time_samples = zeros(nRepeat,1);
        for ir = 1:nRepeat
            [Err_L1(kk),time_samples(ir)] = CPU_time_cost_L1_CLG(...
                kappa,a,b,Alpha,Beta,Gamma,Ntilde_L1(kk),M_now,...
                current_time,grade);
        end
        CPU_L1(kk) = median(time_samples);
    end

    fprintf('CIM-CLG: N = %g, L2 error = %.6e, CPU = %.6e s\n',...
        N_CIM(kk),Err_CIM(kk),CPU_CIM(kk));
    fprintf('L1-CLG : Ntilde = %g, L2 error = %.6e, CPU = %.6e s\n',...
        Ntilde_L1(kk),Err_L1(kk),CPU_L1(kk));
end

%--------------------------------------------------------------
% Print the numerical results.
%--------------------------------------------------------------
fprintf('\n');
fprintf('==========================================================================\n');
fprintf(' M       N(CIM)       Ntilde(L1)       Error(CIM)       Error(L1)       CPU(CIM)       CPU(L1)\n');
fprintf('==========================================================================\n');
for kk = 1:length(M)
    fprintf('%2d      %6g         %6g        %11.4e    %11.4e    %11.4e    %11.4e\n',...
        M(kk),N_CIM(kk),Ntilde_L1(kk),Err_CIM(kk),Err_L1(kk),...
        CPU_CIM(kk),CPU_L1(kk));
end
fprintf('==========================================================================\n');

%--------------------------------------------------------------
% Figure 1: temporal discretization counts versus M.
%--------------------------------------------------------------
figure(1); clf;
set(gcf,'Position',[100,100,720,520],'Color','w');
plot(M,N_CIM,'o-','LineWidth',1.6,'MarkerSize',7,...
     'MarkerFaceColor',[0.12,0.47,0.71]);
hold on;
plot(M,Ntilde_L1,'s-','LineWidth',1.6,'MarkerSize',7,...
     'MarkerFaceColor',[0.85,0.33,0.10]);
grid on; box on;
xlabel('Polynomial degree $M$','Interpreter','latex','FontSize',14);
ylabel('Temporal discretization count','FontSize',14);
title('1D manufactured solution: temporal discretization counts',...
      'FontSize',14);
legend('CIM contour nodes $N$','L1 time steps $\widetilde{N}$',...
       'Interpreter','latex','Location','northwest');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);

%--------------------------------------------------------------
% Figure 2: CPU time versus M.
%--------------------------------------------------------------
figure(2); clf;
set(gcf,'Position',[140,140,720,520],'Color','w');
plot(M,CPU_CIM,'o-','LineWidth',1.6,'MarkerSize',7,...
     'MarkerFaceColor',[0.12,0.47,0.71]);
hold on;
plot(M,CPU_L1,'s-','LineWidth',1.6,'MarkerSize',7,...
     'MarkerFaceColor',[0.85,0.33,0.10]);
grid on; box on;
xlabel('Polynomial degree $M$','Interpreter','latex','FontSize',14);
ylabel('CPU time (s)','FontSize',14);
title(sprintf('1D manufactured solution: CPU time (Tol = %.4e)',tol),...
      'FontSize',14);
legend('CIM-CLG','L1-CLG','Location','northwest');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);

fprintf('\nTwo MATLAB figures have been generated.\n');
fprintf('Figure 1: M versus N and Ntilde.\n');
fprintf('Figure 2: M versus CPU time.\n');
fprintf('No data files are saved automatically.\n');
