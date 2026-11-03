function fig=J9_plot_fft_comparison(R,measurement)
% J9_PLOT_FFT_COMPARISON  Journal-style standalone FFT/direct comparison.
% Calls produce ONE figure each, never subplots and never save files.
%
% 1D: polynomial degrees [16 32 64 128 256], actual Example 2.
% 2D: polynomial degrees [8 12 16 20 24], actual six-peak K=J=2.
% The colors, font sizes, bar widths, and figure dimensions are identical.
%
% Use 'total' for the complete solver CPU/wall-clock comparison (default).
% Use 'transform' for the isolated DCT-I cost (optional).

if nargin<2, measurement='total'; end
measurement=lower(char(measurement));

switch measurement
    case 'total'
        Y=R.total;
        suffix='total runtime';
        ylabelText='Total elapsed time (s)';
    case 'transform'
        Y=R.dct;
        suffix='CGL transform runtime';
        ylabelText='DCT-I elapsed time (s)';
    otherwise
        error('measurement must be total or transform');
end

if R.dimension==1
    caption=sprintf('1D CIM--CLG: %s ($N=%d$, $t=%.2f$)', ...
        suffix,R.N,R.t);
    windowTitle='Comment 9: 1D Example 2 FFT vs direct';
elseif R.dimension==2
    caption=sprintf('2D CIM--CLG: %s ($K=J=2$, $N=%d$)', ...
        suffix,R.N);
    windowTitle='Comment 9: 2D six-peak FFT vs direct';
else
    error('Unknown spatial dimension');
end

fig=figure('Color','w','Units','pixels', ...
    'Position',[150 150 725 555], ...
    'Name',windowTitle,'NumberTitle','off');

ax=axes('Parent',fig,'Units','normalized', ...
    'Position',[0.092 0.135 0.872 0.788], ...
    'FontName','Times New Roman', ...
    'FontSize',13,'LineWidth',1.0, ...
    'Box','on','TickDir','out');
x=(1:numel(R.M)).';

bars=bar(ax,x,Y,0.68,'grouped');
bars(1).FaceColor=[31 119 180]/255;  % blue: FFT
bars(2).FaceColor=[255 127 14]/255;  % orange: no FFT
bars(1).EdgeColor='none';
bars(2).EdgeColor='none';

ax.XLim=[0.48,numel(R.M)+0.52];
ax.XTick=x;
ax.XTickLabel=arrayfun(@num2str,R.M(:),'UniformOutput',false);
ax.XGrid='off';
ax.YGrid='on';
ax.GridAlpha=0.15;

xlabel(ax,'Polynomial degree $M$', ...
    'Interpreter','latex','FontSize',18);
ylabel(ax,ylabelText,'Interpreter','latex','FontSize',18);
title(ax,caption, ...
    'Interpreter','latex', ...
    'FontSize',17,'FontWeight','normal');

lg=legend(ax,bars,{ ...
    'FFT-based CGL (DCT-I)', ...
    'Direct CGL (DCT-I; no FFT)'}, ...
    'Interpreter','latex', ...
    'Location','northwest', ...
    'FontSize',12,'Box','on');
lg.Color='white';
lg.EdgeColor=[0.72 0.72 0.72];
if isprop(lg,'ItemTokenSize')
    lg.ItemTokenSize=[30 14];
end

maximum=max(Y(:));
if maximum>0
    ylim(ax,[0,1.12*maximum]);
end
drawnow;
end
