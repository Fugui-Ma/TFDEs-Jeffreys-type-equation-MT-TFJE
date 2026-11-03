function R=Main_Jeffreys_compare_2d(showTransform)
% 2D Example 3: compare FULL CIM--CLG runtimes using FFT/direct CGL DCT-I.
%
% Run:
%   Jeffreys_compare_2d             % ONE independent total-time figure
%   Jeffreys_compare_2d(true)       % plus a separate DCT-only figure
%
% The initial condition is the six-peak Dirichlet-compatible field:
% p0(x,y)=(1-x^2)(1-y^2)*sum_{k=0}^5
%         exp(-64*((x-.8*cos(k*pi/3))^2+(y-.8*sin(k*pi/3))^2)).
% This is HOMOGENEOUS f=0, K=J=2, so initial data are transformed ONCE.
% All N contour systems of size (M-1)^2 are solved by A\rhs.
% NO data files or figures are saved.
if nargin<1, showTransform=false; end

cfg.problem='2D';
cfg.dimension=2;
cfg.caseName='sixpeak';
% Note: 2D has (M-1)^2 DOFs; using M=256 would be excessively expensive.
% The largest case below M=24 has 529 spatial unknowns.
cfg.Mlist=[16 32 64 128 256];
cfg.N=50;
cfg.t=0.50;
cfg.repeats=5;

cfg.p.alpha=0.50;
cfg.p.beta=0.35;
cfg.p.gamma=0.45;
cfg.p.a=1.0;
cfg.p.b=0.5;
cfg.p.ax=[0.10 0.18];
cfg.p.aa=[0.35 0.20];
cfg.p.bx=[0.05 0.10];
cfg.p.bb=[0.25 0.15];
cfg.p.kappa=2;   % not used when f=0 with six-peak initial data
cfg.p.t0=0.25;
cfg.p.Lambda=3.0;

R=J9_fft_direct_benchmark(cfg);
J9_plot_fft_comparison(R,'total');
if showTransform
    J9_plot_fft_comparison(R,'transform');
end
end
