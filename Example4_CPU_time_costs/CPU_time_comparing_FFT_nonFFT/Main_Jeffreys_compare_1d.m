function R=Jeffreys_compare_1d(showTransform)
%**************************************************************************
%**    Program： By using the CIM-CGL method to solve the 2D Case :
%                with zero Dirichlet bundary condition
%            
%Explainations:
%   this expriment results showing in  Fig7
%***Author：
%           F.G. MA
%*****Finishi Data：
%           2026-10-23
%**************************************************************************
% 1D Example 2: compare FULL CIM--CLG runtimes using FFT/direct CGL DCT-I.
%
% Run:
%   Jeffreys_compare_1d             % ONE independent total-time figure
%   Jeffreys_compare_1d(true)       % plus a separate DCT-only figure
%
% Exact solution:
% p(t,x)=t^kappa (1.1+sin(pi*x))(1-x^2), with K=J=0.
% The full solver evaluates fhat at EVERY contour node and uses A\rhs.
% NO data files or figures are saved.
if nargin<1, showTransform=false; end

cfg.problem='1D';
cfg.dimension=1;
cfg.caseName='analytic';
cfg.Mlist=[16 32 64 128 256];
cfg.N=50;
cfg.t=0.50;
cfg.repeats=5;

cfg.p.alpha=0.50;
cfg.p.beta=0.35;
cfg.p.gamma=0.45;
cfg.p.a=1.0;
cfg.p.b=0.5;
cfg.p.ax=[]; cfg.p.aa=[];
cfg.p.bx=[]; cfg.p.bb=[];
cfg.p.kappa=0.25;
cfg.p.t0=0.25;
cfg.p.Lambda=3.0;

R=J9_fft_direct_benchmark(cfg);
J9_plot_fft_comparison(R,'total');
if showTransform
    J9_plot_fft_comparison(R,'transform');
end
end
