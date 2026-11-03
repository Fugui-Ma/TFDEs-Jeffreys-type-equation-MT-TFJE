function [fig2D,fig3D] = Plot_Jeffreys2D(nodes,U,t,Alpha,Beta,Gamma)
%PLOT_JEFFREYS2D Clear, bright, publication-style visualization.
%   Plot_Jeffreys2D(nodes,U,t,Alpha,Beta,Gamma)
%
% nodes: ascending 1D LGL nodes in [-1,1].
% U:     matrix on the meshgrid(nodes,nodes), rows are y, columns are x.
%
% This function ONLY draws the supplied numerical solution. It does not
% modify the Jeffreys equation, Galerkin coefficients, or CIM approximation.
%
% MATLAB R2025a. Uses built-in functionality only.

    nodes = nodes(:);
    U = real(double(U));
    if numel(nodes) ~= size(U,1) || size(U,1) ~= size(U,2)
        error('size(U) must equal [numel(nodes), numel(nodes)].');
    end
    if any(diff(nodes)<=0)
        error('nodes must be strictly increasing.');
    end
    if any(~isfinite(U(:)))
        error('U contains NaN or Inf; check the numerical solver.');
    end

    % Smooth interpolated DISPLAY grid. U itself remains unchanged.
    fine = linspace(nodes(1),nodes(end),401);
    [xx,yy] = meshgrid(fine,fine);
    Z = interp2(nodes,nodes,U,xx,yy,'makima');

    zmin=min(Z(:));
    zmax=max(Z(:));
    if zmax == zmin
        zmin=zmin-0.5; zmax=zmax+0.5;
    end

    % Keep true data limits. In particular, do not hide negative numerical
    % values by replacing the computed solution with max(U,0).
    lower=min(0,zmin);
    upper=zmax;
    if upper <= lower
        upper=lower+eps(max(1,abs(lower)));
    end

    % Bright, balanced palette: indigo -> blue -> cyan -> lemon -> orange -> red.
    anchors = [ ...
        0.08 0.10 0.42; ...
        0.16 0.36 0.89; ...
        0.09 0.76 0.90; ...
        0.64 0.93 0.64; ...
        1.00 0.91 0.34; ...
        1.00 0.54 0.17; ...
        0.84 0.13 0.28];
    cmap = interp1(linspace(0,1,size(anchors,1)), ...
        anchors,linspace(0,1,256),'pchip');
    cmap = max(0,min(1,cmap));

    titletext = sprintf( ...
        '$\\alpha=%.2f,\\quad\\beta=%.2f,\\quad\\gamma=%.2f,\\quad t=%g$', ...
        Alpha,Beta,Gamma,t);

    % ------------------ Contour plot ------------------
    fig2D=figure('Color','w','Units','pixels','Position',[80 90 1020 750]);
    ax=axes('Parent',fig2D,'Position',[.10 .11 .73 .79]);

    levels=linspace(lower,upper,46);
    contourf(ax,xx,yy,Z,levels,'LineStyle','none');
    hold(ax,'on');

    % Sparse subtle isolines, instead of 20 heavy black contours.
    thinlevels=linspace(lower,upper,9);
    contour(ax,xx,yy,Z,thinlevels(2:end-1), ...
            'LineColor',[0.12 0.22 0.33],'LineWidth',0.35);

    colormap(ax,cmap);
    clim(ax,[lower upper]);
    axis(ax,'equal');
    xlim(ax,[nodes(1),nodes(end)]);
    ylim(ax,[nodes(1),nodes(end)]);
    set(ax,'FontName','Times New Roman','FontSize',14,...
        'LineWidth',1.15,'TickDir','out','Box','on','Layer','top');
    ax.XTick=linspace(-1,1,5);
    ax.YTick=linspace(-1,1,5);

    xlabel(ax,'$x$','Interpreter','latex','FontSize',19);
    ylabel(ax,'$y$','Interpreter','latex','FontSize',19);
    title(ax,titletext,'Interpreter','latex',...
        'FontSize',18,'FontWeight','normal');

    cb=colorbar(ax,'eastoutside');
    cb.Position=[.86 .17 .025 .69];
    cb.FontName='Times New Roman';
    cb.FontSize=12;
    cb.LineWidth=0.8;
    cb.TickDirection='out';
    cb.Ticks=linspace(lower,upper,6);
    cb.Label.String='$p(t,x,y)$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=16;

    % ------------------ 3D surface --------------------
    fig3D=figure('Color','w','Units','pixels','Position',[1140 90 930 735]);
    ax3=axes('Parent',fig3D,'Position',[.08 .12 .75 .77]);
    surf(ax3,xx,yy,Z,'EdgeColor','none','FaceColor','interp');
    shading(ax3,'interp');
    colormap(ax3,cmap);
    clim(ax3,[lower upper]);
    view(ax3,-53,33);
    xlim(ax3,[-1,1]);
    ylim(ax3,[-1,1]);
    axis(ax3,'tight');
    set(ax3,'FontName','Times New Roman','FontSize',13,...
        'LineWidth',1.0,'TickDir','out','Box','on');
    xlabel(ax3,'$x$','Interpreter','latex','FontSize',17);
    ylabel(ax3,'$y$','Interpreter','latex','FontSize',17);
    zlabel(ax3,'$p(t,x,y)$','Interpreter','latex','FontSize',17);
    title(ax3,titletext,'Interpreter','latex','FontSize',17);
    c3=colorbar(ax3);
    c3.FontSize=12;
    c3.Label.String='$p(t,x,y)$';
    c3.Label.Interpreter='latex';

    % Un-comment to export publication-ready figures as PNG and PDF:
    % exportgraphics(fig2D,'Jeffreys_2D_Contour.png','Resolution',400);
    % exportgraphics(fig2D,'Jeffreys_2D_Contour.pdf','ContentType','vector');
    % exportgraphics(fig3D,'Jeffreys_2D_Surface.png','Resolution',400);
end
