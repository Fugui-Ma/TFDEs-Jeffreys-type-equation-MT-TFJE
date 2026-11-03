function R=J9_fft_direct_benchmark(cfg)
% C9_FFT_DIRECT_BENCHMARK  Fair end-to-end FFT vs direct CGL benchmark.
%
% At each polynomial degree M:
%   1. Warm up both methods.
%   2. Alternate their order on successive repetitions to limit order bias.
%   3. Select each method's median full wall-clock-time run.
%   4. Compare the final Galerkin coefficients for numerical equivalence.
%
% The methods use IDENTICAL model, contour, Galerkin matrices and A\rhs
% solver. Only FFT-based versus direct DCT-I changes.
% The c9_source_solver.m subroutine is the solver taken from the user's
% latest Comment9 code package; it returns measured, not invented, times.
%
% R.total(:,1)=FFT; R.total(:,2)=Direct cosine transform.
% R.dct(:,1)=FFT DCT-I time; R.dct(:,2)=Direct DCT-I time.
% No file output; timing results are only returned in R/printed.

Mlist=cfg.Mlist(:).';
n=numel(Mlist);
methods={'fft','direct'};

R.problem=cfg.problem;
R.M=Mlist;
R.N=cfg.N;
R.t=cfg.t;
R.dimension=cfg.dimension;
R.caseName=cfg.caseName;
R.repeats=cfg.repeats;
R.total=zeros(n,2);
R.dct=zeros(n,2);
R.stages=zeros(n,2,8);
R.relativeDifference=zeros(n,1);
% speedup >1 means FFT is faster; speedup <1 means direct is faster.
R.speedup=zeros(n,1);
R.L2error=nan(n,2);
R.dof=zeros(n,1);

assert(cfg.repeats>=3 && cfg.repeats==round(cfg.repeats), ...
    'Use three or more integer timing repetitions.');
assert(cfg.N>=2 && cfg.N==round(cfg.N),'N must be an integer >=2.');

fprintf('\n%s: FFT versus direct CGL; N=%d; t=%.3f; Q_t=1\n', ...
    cfg.problem,cfg.N,cfg.t);
fprintf('All full runtimes include setup, source/RHS, matrix assembly, and N solves by A\\b.\n');

for i=1:n
    M=Mlist(i);
    dof=(M-1)^cfg.dimension;
    R.dof(i)=dof;
    p=cfg.p;

    % MATLAB JIT / FFT planning warm-up, excluded from timing.
    for methodId=1:2
        J9_source_solver(M,cfg.N,cfg.t,cfg.dimension, ...
            cfg.caseName,p,methods{methodId});
    end

    % Store one full solver output per timed run, for each transform.
    attempts=cell(cfg.repeats,2);
    wall=zeros(cfg.repeats,2);
    for rep=1:cfg.repeats
        if mod(rep,2)==1
            order=[1 2]; % odd runs: FFT then direct
        else
            order=[2 1]; % even runs: direct then FFT
        end
        for j=order
            run=J9_source_solver(M,cfg.N,cfg.t,cfg.dimension, ...
                cfg.caseName,p,methods{j});
            attempts{rep,j}=run;
            wall(rep,j)=run.T.total;
        end
    end

    representative=cell(1,2);
    for j=1:2
        target=median(wall(:,j));
        [~,idx]=min(abs(wall(:,j)-target));
        representative{j}=attempts{idx,j};
        R.total(i,j)=representative{j}.T.total;
        R.dct(i,j)=representative{j}.T.dct;
        R.stages(i,j,:)=reshape(representative{j}.T.stage,1,1,8);

        % Verify that the FULL solver timing includes every stage.
        assert(abs(sum(representative{j}.T.stage) ...
               -representative{j}.T.total)< ...
               1e-8+1e-8*representative{j}.T.total, ...
            'Eight phase timers do not reproduce the full elapsed time.');
    end

    R.speedup(i)=R.total(i,2)/max(R.total(i,1),eps);

    uFFT=representative{1}.coeff;
    uDirect=representative{2}.coeff;
    R.relativeDifference(i)=norm(uFFT(:)-uDirect(:)) ...
        /max(norm(uFFT(:)),eps);

    if R.relativeDifference(i)>1e-8
        error(['Numerical mismatch at M=%d: relative difference %.3e. ',...
               'Check the FFT/direct DCT-I implementations.'], ...
              M,R.relativeDifference(i));
    end

    if cfg.dimension==1
        % Error to the exact solution is NOT inside the timed interval.
        R.L2error(i,1)=J9_source_error(uFFT,M,1,'analytic',p,cfg.t);
        R.L2error(i,2)=J9_source_error(uDirect,M,1,'analytic',p,cfg.t);
        fprintf('M=%3d | dof=%4d | total FFT %.5f s, direct %.5f s | ',...
            M,dof,R.total(i,1),R.total(i,2));
        fprintf('DCT FFT %.5f s, direct %.5f s | L2 FFT %.2e | diff %.2e\n', ...
            R.dct(i,1),R.dct(i,2), ...
            R.L2error(i,1),R.relativeDifference(i));
        fprintf('         direct/FFT full runtime ratio: %.3f\n',R.speedup(i));
    else
        fprintf('M=%3d | dof=%4d | total FFT %.5f s, direct %.5f s | ', ...
            M,dof,R.total(i,1),R.total(i,2));
        fprintf('initial DCT FFT %.5f s, direct %.5f s | diff %.2e\n', ...
            R.dct(i,1),R.dct(i,2),R.relativeDifference(i));
        fprintf('         direct/FFT full runtime ratio: %.3f\n',R.speedup(i));
    end
end

if cfg.dimension==2
    fprintf(['Note: f=0 in 2D. Initial-data DCT-I is computed only ONCE, ',...
             'so FFT need not substantially reduce total 2D runtime.\n']);
end
fprintf('Storage per output coefficient vector (real double): 8*(M-1)^d bytes.\n');
fprintf('Streaming avoids storing N complex contour coefficient vectors.\n\n');
end
