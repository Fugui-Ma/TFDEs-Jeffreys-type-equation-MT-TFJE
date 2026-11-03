function J9_source_plot(kind,varargin)
% C9_SOURCE_PLOT: one figure per physical problem, never saves files.
%
% Timing columns in the order returned by c9_source_solver:
%   1 setup and spatial data; 
%   2 fhat/symbol evaluation;
%   3 CGL transform; 
%   4 Galerkin RHS vector assembly;
%   5 Galerkin matrix assembly; 
%   6 complete linear solve (A\b);
%   7 CIM accumulation; 
%   8 solution reconstruction.
%
% Usage:
%   c9_source_plot('breakdown',N,stages,'1d')
%   c9_source_plot('breakdown',N,stages,'2d')
%   c9_source_plot('paired',N,FFTstages,DirectStages,'1d'/'2d')
%   c9_source_plot('compare',M,totalTimes,dctOnlyTimes)
%   c9_source_plot('spatial',M,errors)

switch lower(char(kind))
    case 'breakdown'
        N=varargin{1};
        T=varargin{2};
        problem=lower(char(varargin{3}));
        assert(size(T,2)==8,'Eight timing phases are required.');
        [fig,ax]=newFigure(problem,'CIM--CLG runtime breakdown');
        drawStack(ax,N,T,problem,false);
        set(fig,'Color','w');

    case 'paired'
        N=varargin{1};F=varargin{2};D=varargin{3};
        problem=lower(char(varargin{4}));
        assert(isequal(size(F),size(D)) && size(F,2)==8);
        [~,ax]=newFigure(problem,'FFT versus direct CGL transforms');
        drawPaired(ax,N,F,D,problem);

    case 'compare'
        M=varargin{1};whole=varargin{2};dct=varargin{3};
        drawGrouped(M,whole,'1D Example 2: total CIM--CLG runtime (N=50)', ...
            'Total CPU time (s)');
        drawGrouped(M,dct,'1D Example 2: CGL DCT-I runtime (N=50)', ...
            'CGL transform CPU time (s)');

    case 'spatial'
        M=varargin{1};E=varargin{2};
        fig=figure('Color','w','Units','pixels', ...
            'Position',[150 150 820 555], ...
            'Name','Comment 9: spatial convergence', ...
            'NumberTitle','off');
        ax=axes('Parent',fig,'Position',[.105 .17 .855 .745]);
        semilogy(ax,M,E(:,1),'-o','LineWidth',1.7,'MarkerSize',7);hold(ax,'on');
        semilogy(ax,M,E(:,2),'-s','LineWidth',1.7,'MarkerSize',7);hold(ax,'off');
        xlabel(ax,'Polynomial degree M');
        ylabel(ax,'L^2 spatial error');
        title(ax,'Spatial convergence: analytic versus finite regularity');
        legend(ax,'Analytic spatial factor','Finite spatial regularity', ...
            'Location','southwest');
        grid(ax,'on');
        set(ax,'FontName','Arial','FontSize',13,'Box','on');

    otherwise
        error('Unknown c9_source_plot mode: %s',char(kind));
end
end

function [fig,ax]=newFigure(problem,purpose)
% Identical sizes for 1D and 2D. Recreates the proportions of the
% user's preferred original panel, now as a standalone figure.
if strcmpi(problem,'1d')
    name='1D Example 2';
elseif strcmpi(problem,'2d')
    name='2D Example 3 (K=J=2)';
else
    error('problem must be 1d or 2d.');
end
fig=figure('Color','w','Units','pixels', ...
    'Position',[150 150 725 555], ...
    'Name',[name,' - ',purpose], ...
    'NumberTitle','off');
ax=axes('Parent',fig,'Units','normalized', ...
    'Position',[.105 .155 .865 .780], ...
    'FontName','Arial','FontSize',18,'Box','on', ...
    'LineWidth',1.0,'TickDir','out');
end

function colors=phaseColors()
% Retains the original blue/orange/green/purple accents while distinguishing
% the extra setup, transformation, projection and reconstruction work.
colors=[ ...
     99 110 118;  % 1 setup: gray
     31 119 180;  % 2 transformed source or symbols: blue
     23 190 207;  % 3 DCT-I: cyan
    188 189  34;  % 4 source/RHS vector assembly: olive
    255 127  14;  % 5 Galerkin matrix assembly: orange
     44 160  44;  % 6 MATLAB A\b linear solve: green
    148 103 189;  % 7 CIM accumulation: purple
    140  86  75] / 255; % 8 physical reconstruction: brown
end

function names = stageNames(problem)

if strcmpi(problem,'1d')
    names = {
        'Basis and contour setup', ...
        'Source evaluation: $\widehat{f}(z_k)$', ...
        'CGL transform of RHS', ...
        'RHS vector assembly: $\mathbf{b}_k$', ...
        'Galerkin matrix assembly: $\mathbf{A}(z_k)$', ...
        'Linear system solve: $\mathbf{A}(z_k)\mathbf{c}_k=\mathbf{b}_k$', ...
        'CIM quadrature accumulation', ...
        'Solution reconstruction'
    };
else
    names = {
        'Basis and initial-data setup', ...
        'Laplace-domain initial RHS', ...
        'CGL transform of $p_0$', ...
        'Initial RHS assembly: $\mathbf{b}_{0,M}$', ...
        'Galerkin matrix assembly: $\mathbf{A}(z_k)$', ...
        'Linear system solve: $\mathbf{A}(z_k)\mathbf{c}_k=\mathbf{b}_k$', ...
        'CIM quadrature accumulation', ...
        'Solution reconstruction'
    };
end

end

function drawStack(ax,N,T,problem,isPaired)
N=N(:).';
assert(size(T,1)==numel(N) && size(T,2)==8);
assert(all(isfinite(T(:))) && all(T(:)>=-1e-12), ...
    'Timing data are not valid.');

c=phaseColors();
h=bar(ax,N,max(T,0),.15,'stacked');
for k=1:8
    h(k).FaceColor=c(k,:);
    h(k).EdgeColor='none';
end
formatAxes(ax,N,problem,isPaired);
makeLegend(ax,h,stageNames(problem));

maxTime=max(sum(T,2));
ylim(ax,[0,max(1e-8,1.09*maxTime)]);
end

function drawPaired(ax,N,F,D,problem)
N=N(:);
assert(size(F,1)==numel(N) && size(F,2)==8);
assert(isequal(size(F),size(D)));
x=zeros(2*numel(N),1);
T=zeros(2*numel(N),8);
for k=1:numel(N)
    x(2*k-1)=N(k)-1.50;    % FFT on the left
    x(2*k)=N(k)+1.50;      % direct cosine transform on the right
    T(2*k-1,:)=F(k,:);
    T(2*k,:)=D(k,:);
end
c=phaseColors();
h=bar(ax,x,max(T,0),.15,'stacked');
for k=1:8
    h(k).FaceColor=c(k,:);
    h(k).EdgeColor='none';
end
hold(ax,'on');
for k=1:numel(N)
    rectangle('Parent',ax, ...
        'Position',[x(2*k)-.97,0,1.94,sum(D(k,:))], ...
        'FaceColor','none', ...
        'EdgeColor',[.16 .16 .16],'LineWidth',.8);
end
hold(ax,'off');
formatAxes(ax,N,problem,true);
makeLegend(ax,h,stageNames(problem));
ylim(ax,[0,max(1e-8,1.09*max(sum(T,2)))]);
end

function formatAxes(ax,N,problem,isPaired)
if strcmpi(problem,'1d')
    head='1D CIM-CLG ($M=24$)';
else
    head='2D CIM-CLG ($M=24)$';
end
if isPaired
    head=[head,' - FFT left, direct right'];
end
set(ax,'XLim',[8 44],'XTick',10:15:40,'Layer','top', ...
    'FontName','Arial','FontSize',18,'Box','on');
xlabel(ax,'Contour nodes $N$', 'Interpreter','latex','FontSize',18);
ylabel(ax,'CPU time (s)', 'Interpreter','latex','FontSize',18);
title(ax,head,'FontSize',16,'FontWeight','normal', ...
    'Interpreter','latex');
grid(ax,'off');
end

function makeLegend(ax,handles,names)
% Compact two-column legend. Clear numerical-stage terminology.
lg = legend(ax,handles,names,'Location','northwest', ...
    'Interpreter','latex', ...
    'FontSize',12, ...
    'Box','on');
lg.FontName='Arial';
lg.Color='white';
lg.EdgeColor=[.80 .80 .80];
if isprop(lg,'NumColumns'),lg.NumColumns=2;end
if isprop(lg,'ItemTokenSize'),lg.ItemTokenSize=[22 11];end
end

function drawGrouped(M,T,head,yText)
figure('Color','w','Units','pixels', ...
    'Position',[100 100 725 655], ...
    'Name',head,'NumberTitle','off');
ax=axes('Position',[.105 .155 .865 .780]);
x=1:numel(M);
b=bar(ax,x,T,.15,'grouped');
b(1).FaceColor=[31 119 180]/255;
b(2).FaceColor=[255 127 14]/255;
b(1).EdgeColor='none';
b(2).EdgeColor='none';
set(ax,'XTick',x,'XTickLabel', ...
    arrayfun(@num2str,M(:).','UniformOutput',false), ...
    'FontName','Arial','FontSize',16,'Box','on');
xlabel(ax,'Polynomial degree M','FontSize',18);
ylabel(ax,yText,'FontSize',18);
title(ax,head,'FontWeight','normal','FontSize',16);
ylim(ax,[0,1.12*max(T(:))+eps]);
legend(ax,b,{'FFT-based DCT-I','Direct cosine transform'}, ...
    'Location','northwest','FontSize',14);
grid(ax,'off');
end
