function Run_Jeffreys2D_KJ2_ST_order()
%**************************************************************************
%**    Program： By using the CIM-CGL method to solve the 2D Case :
%                with zero Dirichlet bundary condition
%            
%Explainations:
%   this expriment results showing in  Table 4
%***Author：
%           F.G. MA
%*****Finishi Data：
%           2026-10-22
%**************************************************************************
% MATLAB Command Window ONLY: K=J=2 CIM--CLG errors and observed orders.
% No tables, graphics, .mat, .csv, .tex, or any files are created.
% Uses the original uploaded six-peak initial function without modification.
%
% Time error: max over sampled times of ||p_M^N - p_M^Nref||_L2.
% Space error: ||p_M^N - p_2M^N||_L2 / sampled L-infinity.
% Orders: log(E_previous/E_current)/log(level_current/level_previous).
% They are empirical reference/doubling orders, NOT exact-solution errors.
clear all; close all; clc

%% ------------------- Parameters (keep user settings) -------------------
N_list = [4 6 10 20 30 50];
N_reference = 100;
M_temporal = 12;
M_space = [6 12 24];
N_spatial = 50;
current_time = 0.5;

% These match the user's last uploaded main program (Lambda=500).
% For the manuscript's Lambda=150 setting, change both values to 150.
t0_time = 0.01;     B_gamma_time = 500;
t0_space = 0.01;    B_gamma_space = 500;

% Preserve the existing six-Gaussian initial function, including pi^2.
initial_function = @Initial_values_sixpeak_KJ2;
num_time_samples = 121; % Sampled time-interval maximum, not exact supremum.
num_linf_grid = 401;    % Sampled max norm for the spatial error.

triples = [0.25 0.55 0.75;
           0.50 0.35 0.45;
           0.75 0.15 0.15];
case_names = {'a','b','c'};

%% --------------------------- Input checks ------------------------------
assert(all(diff(N_list)>0) && all(N_list>0) && ...
       all(N_list==round(N_list)),'N_list must be increasing integers.');
assert(N_reference>max(N_list),'N_reference must exceed all N_list.');
assert(all(diff(M_space)>0) && all(M_space>=4) && ...
       all(M_space==round(M_space)),'M_space must be increasing integers.');
assert(current_time>=t0_time && current_time<=t0_time*B_gamma_time, ...
       'Observation time must lie in the temporal CIM window.');
assert(current_time>=t0_space && current_time<=t0_space*B_gamma_space, ...
       'Observation time must lie in the spatial CIM window.');
assert(num_time_samples>=21 && num_linf_grid>=101);
assert(size(triples,1)==3);

base = struct('K',2,'J',2,'a',1,'b',0.5, ...
    'alpha_terms',[0.10,0.18], 'beta_terms',[0.05,0.10], ...
    'a_terms',[0.35,0.20], 'b_terms',[0.25,0.15], ...
    'N',N_spatial,'t',current_time, ...
    't0',t0_space,'B_gamma',B_gamma_space);

%% -------------------- Assemble spatial systems -------------------------
spaceDegrees = unique([M_space,2*M_space]);
spatial = cell(1,numel(spaceDegrees));

fprintf('\n===============================================================\n');
fprintf('2D Jeffreys CIM--CLG: K=J=2, Command Window results ONLY\n');
fprintf('Initial function: %s\n',func2str(initial_function));
fprintf('a=%.3g, b=%.3g, t*=%.3g; alpha_k=[.10 .18]; beta_j=[.05 .10]\n', ...
    base.a,base.b,current_time);
fprintf('a_k=[.35 .20]; b_j=[.25 .15]\n');
fprintf('Time: M=%d, Nref=%d, t0=%.3g, Lambda=%.3g\n', ...
    M_temporal,N_reference,t0_time,B_gamma_time);
fprintf('Space: N=%d, M=%s, t0=%.3g, Lambda=%.3g\n', ...
    N_spatial,mat2str(M_space),t0_space,B_gamma_space);
fprintf('No data or figure files will be written.\n');
fprintf('===============================================================\n');

for im=1:numel(spaceDegrees)
    M=spaceDegrees(im);
    fprintf('Building Galerkin system: M=%d\n',M);
    spatial{im}=Build_Jeffreys2D_Spatial(M,initial_function);
end
idxT=find(spaceDegrees==M_temporal,1);
if isempty(idxT)
    spT=Build_Jeffreys2D_Spatial(M_temporal,initial_function);
else
    spT=spatial{idxT};
end

% The manuscript defines error_N using a maximum over [t0,Lambda*t0].
% Include t* exactly for the separate fixed-time comparison.
times=unique([logspace(log10(t0_time), ...
             log10(t0_time*B_gamma_time),num_time_samples),current_time]);
[~,obsIdx]=min(abs(times-current_time));
assert(abs(times(obsIdx)-current_time) < 1e-12, ...
       'Observation time was not included on the time grid.');
denseGrid=linspace(-1,1,num_linf_grid).';

%% ------------------- Solve and print three cases -----------------------
for ic=1:3
    cfg=base;
    cfg.Alpha=triples(ic,1);
    cfg.Beta=triples(ic,2);
    cfg.Gamma=triples(ic,3);

    fprintf('\n---------------------------------------------------------------\n');
    fprintf('CASE (%s): alpha=%.2f, beta=%.2f, gamma=%.2f\n', ...
        case_names{ic},cfg.Alpha,cfg.Beta,cfg.Gamma);
    fprintf('---------------------------------------------------------------\n');

    % [1] TIME ERRORS, same Galerkin space and multiple contour points N.
    cfgTime=cfg;
    cfgTime.t0=t0_time;
    cfgTime.B_gamma=B_gamma_time;
    referenceModes=Jeffreys2D_CIM_ModalTimes_KJ2( ...
        spT,cfgTime,N_reference,times);

    errorMax=zeros(numel(N_list),1);
    errorAtT=zeros(numel(N_list),1);
    for jn=1:numel(N_list)
        modesN=Jeffreys2D_CIM_ModalTimes_KJ2( ...
            spT,cfgTime,N_list(jn),times);
        differences=referenceModes-modesN;
        errs=sqrt(sum(abs(differences).^2,1));
        errorMax(jn)=max(errs);
        errorAtT(jn)=errs(obsIdx);
    end

    % Time convergence Order_t as defined in the user's manuscript.
    [orderMax,~]=Jeffreys2D_ObservedOrders_KJ2(N_list,errorMax);
    [orderAtT,~]=Jeffreys2D_ObservedOrders_KJ2(N_list,errorAtT);

    fprintf('\nTIME convergence (M=%d, Nref=%d):\n',M_temporal,N_reference);
    fprintf('  error_N(max): sampled maximum on [%.3g, %.3g], %d points\n', ...
        t0_time,t0_time*B_gamma_time,numel(times));
    fprintf('  error_N(t*): L2 difference at t*=%.3g\n',current_time);
    for jn=1:numel(N_list)
        fprintf('  N=%3d | error_N(max)=%.5E | Order_t=%8s | ', ...
            N_list(jn),errorMax(jn),displayRate(orderMax(jn)));
        fprintf('error_N(t*)=%.5E | Order_t(t*)=%8s\n', ...
            errorAtT(jn),displayRate(orderAtT(jn)));
    end

    % [2] SPACE ERRORS, same time discretization and doubling M -> 2M.
    cfgSpace=cfg;
    cfgSpace.N=N_spatial;
    cfgSpace.t0=t0_space;
    cfgSpace.B_gamma=B_gamma_space;
    results=cell(1,numel(spaceDegrees));
    for im=1:numel(spaceDegrees)
        results{im}=Solve_Jeffreys2D_KJ2(spatial{im},cfgSpace);
    end

    errorL2=zeros(numel(M_space),1);
    errorLinf=zeros(numel(M_space),1);
    for jm=1:numel(M_space)
        M=M_space(jm);
        idxLow=find(spaceDegrees==M,1);
        idxHigh=find(spaceDegrees==2*M,1);
        e=Jeffreys2D_DoublingError_KJ2( ...
            results{idxLow},results{idxHigh},spatial{idxHigh},denseGrid);
        errorL2(jm)=e.L2;
        errorLinf(jm)=e.Linf_sample;
    end
    [orderL2,~]=Jeffreys2D_ObservedOrders_KJ2(M_space,errorL2);
    [orderLinf,~]=Jeffreys2D_ObservedOrders_KJ2(M_space,errorLinf);

    fprintf('\nSPACE convergence (N=%d, t*=%.3g):\n',N_spatial,current_time);
    fprintf('  error_M: p_(2M)^N minus p_M^N; L_inf uses a %dx%d grid\n', ...
        num_linf_grid,num_linf_grid);
    for jm=1:numel(M_space)
        fprintf('  M=%3d | L2 error=%.5E | Order_s(L2)=%8s | ', ...
            M_space(jm),errorL2(jm),displayRate(orderL2(jm)));
        fprintf('Linf error=%.5E | Order_s(Linf)=%8s\n', ...
            errorLinf(jm),displayRate(orderLinf(jm)));
    end
end

fprintf('\n===============================================================\n');
fprintf('Completed. No MAT/CSV/TEX/TXT/FIG/PNG files were created.\n');
fprintf('--- means first refinement level or error near roundoff.\n');
fprintf('All errors are numerical-reference indicators, not exact PDE errors.\n');
fprintf('===============================================================\n');
end

function str=displayRate(r)
% Blank rate on the first row and when too close to floating-point noise.
if ~isfinite(r)
    str='---';
else
    str=sprintf('%.4f',r);
end
end
