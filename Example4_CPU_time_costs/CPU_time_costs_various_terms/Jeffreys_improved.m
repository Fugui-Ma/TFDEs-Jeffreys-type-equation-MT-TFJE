function Jeffreys_improved(reference)
%**************************************************************************
%**    Program： By using the CIM-CGL method to solve the 2D Case :
%                with zero Dirichlet bundary condition
%            
%Explainations:
%   this expriment results showing in  Fig7-8
%***Author：
%           F.G. MA
%*****Finishi Data：
%           2026-10-23
%**************************************************************************
% COMMENT9_IMPROVED - reproducible Comment 9 numerical checks.
%
% especially Jeffreys_verify.m. The CIM contour, Chebyshev interpolation,
% Legendre-Galerkin operators follow that source; each system is now solved by A\b.
%
% Examples:
%    Jeffreys_improved                % 1D and 2D in TWO SEPARATE figure windows
%    Jeffreys_improved('1d')          % only 1D example and figure
%    Jeffreys_improved('2d')          % only 2D six-peak example and figure
%    Jeffreys_improved('paired')      % two separate FFT/direct stage comparisons
%    Jeffreys_improved('compare')     % FFT vs no FFT: two separate figures
%    Jeffreys_improved('spatial')     % spatial regularity convergence check
%    Jeffreys_improved('outputs')     % output-time/memory scaling in console
%    Jeffreys_improved('check')       % DCT, contour, mldivide self-tests
%
% NO CSV, MAT, PNG, or other files are saved.
% Requires these files in the same MATLAB directory:
%    comment9_improved.m, c9_source_solver.m, c9_source_plot.m,
%    J9_source_error.m

if nargin<1 || isempty(mode)
    mode='reference';
end
mode=lower(char(mode));

% The first case is the one-dimensional manufactured Example 2, K=J=0.
p1=parameters_1d();
% The second is the revised two-dimensional six-peak Example 3, K=J=2.
p2=parameters_2d();

switch mode
    case 'reference'
        % Two independent plots, matching the original narrow-bar style:
        % 1D and 2D each have their own figure with N=12/24/40.
        Nlist=[12 24 40];
        M1=24; M2=24;
        t=.5;
        repeats=5;
        fprintf('1D Example 2: K=J=0, M=%d, Q=1, t=%.2f\n',M1,t);
        r1=benchmark(Nlist,M1,t,1,'analytic',p1,'fft',repeats);
        fprintf('2D Example 3: K=J=2, M=%d, %d unknowns, Q=1, t=%.2f\n', ...
            M2,(M2-1)^2,t);
        r2=benchmark(Nlist,M2,t,2,'sixpeak',p2,'fft',repeats);
        J9_source_plot('breakdown',Nlist,r1.stage,'1d');
        J9_source_plot('breakdown',Nlist,r2.stage,'2d');

    case '1d'
        Nlist=[12 24 40];
        r=benchmark(Nlist,24,.5,1,'analytic',p1,'fft',5);
        J9_source_plot('breakdown',Nlist,r.stage,'1d');

    case '2d'
        Nlist=[12 24 40];
        r=benchmark(Nlist,24,.5,2,'sixpeak',p2,'fft',5);
        J9_source_plot('breakdown',Nlist,r.stage,'2d');

    case 'paired'
        % Optional paired variant of the SAME original reference plot.
        % For each N, FFT goes to the left, direct DCT-I to the right.
        Nlist=[12 24 40];
        t=.5;repeats=5;
        f1=benchmark(Nlist,64,t,1,'analytic',p1,'fft',repeats);
        d1=benchmark(Nlist,64,t,1,'analytic',p1,'direct',repeats);
        f2=benchmark(Nlist,24,t,2,'sixpeak',p2,'fft',repeats);
        d2=benchmark(Nlist,24,t,2,'sixpeak',p2,'direct',repeats);
        for i=1:numel(Nlist)
            rel1=norm(f1.coeff{i}-d1.coeff{i})/max(norm(f1.coeff{i}),eps);
            rel2=norm(f2.coeff{i}-d2.coeff{i})/max(norm(f2.coeff{i}),eps);
            assert(rel1<1e-8 && rel2<1e-8, ...
                'FFT and direct methods disagree at N=%d',Nlist(i));
        end
        J9_source_plot('paired',Nlist,f1.stage,d1.stage,'1d');
        J9_source_plot('paired',Nlist,f2.stage,d2.stage,'2d');

    case 'compare'
        % Same PDE and contour parameters for FFT and direct transforms.
        % Figures are independent, as requested earlier.
        Mlist=[16 32 64 128 256];
        N=50;t=.5;repeats=5;
        total=zeros(numel(Mlist),2);
        transform=zeros(numel(Mlist),2);
        fprintf('1D FFT/direct comparison, N=%d, t=%.2f\n',N,t);
        for i=1:numel(Mlist)
            M=Mlist(i);
            rf=benchmark(N,M,t,1,'analytic',p1,'fft',repeats);
            rd=benchmark(N,M,t,1,'analytic',p1,'direct',repeats);
            total(i,:)=[rf.total(1),rd.total(1)];
            transform(i,:)=[rf.dct(1),rd.dct(1)];
            relative=norm(rf.coeff{1}-rd.coeff{1}) ...
                /max(norm(rf.coeff{1}),eps);
            assert(relative<1e-8, ...
                'FFT and non-FFT solutions differ for M=%d',M);
            fprintf('M=%3d | whole solver FFT/direct: %.4f / %.4f s; ', ...
                M,total(i,1),total(i,2));
            fprintf('DCT-I FFT/direct: %.4f / %.4f s; relative diff %.2e\n', ...
                transform(i,1),transform(i,2),relative);
        end
        J9_source_plot('compare',Mlist,total,transform);

    case 'spatial'
        % Optional: spatial convergence under analytic versus finite
        % smoothness. Both are manufactured solutions, p(t,x)=t^2 v(x).
        p=p2;
        p.kappa=2;
        p.t0=.25;p.Lambda=3;
        Mlist=[6 8 10 12 16 20 24];
        N=42;
        E=zeros(numel(Mlist),2);
        for i=1:numel(Mlist)
            a=J9_source_solver(Mlist(i),N,.5,1,'analytic',p,'fft');
            b=J9_source_solver(Mlist(i),N,.5,1,'finite',p,'fft');
            E(i,:)=[ ...
                J9_source_error(a.coeff,Mlist(i),1,'analytic',p,.5), ...
                J9_source_error(b.coeff,Mlist(i),1,'finite',p,.5)];
            fprintf('M=%2d | analytic L2 %.5e | finite-smoothness L2 %.5e\n', ...
                Mlist(i),E(i,1),E(i,2));
        end
        J9_source_plot('spatial',Mlist,E);

    case 'outputs'
        % No additional figure. Print the Q_t cost and storage information
        % explicitly requested by the referee.
        M=24;N=50;
        fprintf('2D six-peak case: M=%d, unknowns=(M-1)^2=%d, N=%d\n', ...
            M,(M-1)^2,N);
        for Q=[1 4 12]
            if Q==1
                times=.5;
            else
                times=linspace(.25,.75,Q);
            end
            r=benchmark(N,M,times,2,'sixpeak',p2,'fft',3);
            fprintf('Q=%2d | wall %.5f s | CIM+reconstruction %.5f s', ...
                Q,r.total(1),r.stage(1,7)+r.stage(1,8));
            fprintf(' | coefficient storage %.3f MiB\n',8*Q*(M-1)^2/1024^2);
        end
        fprintf('Storing all N complex node vectors adds %.3f MiB.\n', ...
            16*N*(M-1)^2/1024^2);

    case 'check'
        J9_source_solver('check');

    otherwise
        error('Unknown mode %s. Use reference, 1d, 2d, paired, compare, spatial, outputs, or check.', mode);
end
end

function p=parameters_1d()
% One-dimensional manufactured Example 2, no additional fractional terms.
p.alpha=.50; p.beta=.35; p.gamma=.45; p.a=1; p.b=.5;
p.ax=[];p.aa=[];p.bx=[];p.bb=[];
p.kappa=.25; p.t0=.25;p.Lambda=3;
end

function p=parameters_2d()
% Revised Example 3, K=J=2, six-peak initial datum.
p.alpha=.50; p.beta=.35; p.gamma=.45; p.a=1; p.b=.5;
p.ax=[.10 .18];p.aa=[.35 .20];
p.bx=[.05 .10];p.bb=[.25 .15];
p.kappa=2;p.t0=.25;p.Lambda=3;
end

function r=benchmark(Nlist,M,times,dim,caseName,p,method,repeats)
Nlist=Nlist(:).';
n=numel(Nlist);
r.stage=zeros(n,8);
r.total=zeros(n,1);
r.dct=zeros(n,1);
r.errorL2=nan(n,1);
r.coeff=cell(n,1);

% Warm-up reduces first-call/JIT and FFT-plan effects.
J9_source_solver(M,Nlist(1),times,dim,caseName,p,method);

for i=1:n
    samples=zeros(repeats,1);
    outputs=cell(repeats,1);
    for k=1:repeats
        out=J9_source_solver(M,Nlist(i),times,dim,caseName,p,method);
        samples(k)=out.T.total;
        outputs{k}=out;
    end
    % Select the run nearest to the median TOTAL time; preserve stages from
    % that SAME run, instead of summing unrelated phase-wise medians.
    target=median(samples);
    [~,j]=min(abs(samples-target));
    chosen=outputs{j};
    r.total(i)=chosen.T.total;
    r.stage(i,:)=chosen.T.stage;
    r.dct(i)=chosen.T.dct;
    if strcmpi(caseName,'sixpeak')
        r.errorL2(i)=NaN;
    else
        r.errorL2(i)=J9_source_error(chosen.coeff,M,dim,caseName,p,times(1));
    end
    r.coeff{i}=chosen.coeff;
    fprintf('  dim=%d M=%3d N=%2d %6s | total %.5f s | DCT %.5f s', ...
        dim,M,Nlist(i),method,r.total(i),r.dct(i));
    if ~isnan(r.errorL2(i))
        fprintf(' | L2 error %.3e',r.errorL2(i));
    end
    fprintf('\n');
    phaseNames={'setup/spatial data','fhat/symbol evaluation','CGL DCT-I', ...
        'RHS vector assembly','matrix assembly','A\b linear solve', ...
        'CIM accumulation','solution reconstruction'};
    fprintf('    detailed phases (s): ');
    for q=1:numel(phaseNames)
        fprintf('%s=%.3g',phaseNames{q},r.stage(i,q));
        if q<numel(phaseNames),fprintf(' | ');end
    end
    fprintf('\n');
end
end
