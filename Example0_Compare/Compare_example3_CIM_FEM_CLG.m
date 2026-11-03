function plot_example3_CIM_FEM_CLG_onecase_t05()
function Run_Jeffreys2D_KJ2_ST_order()
%**************************************************************************
%**    Program： By using the CIM-CGLand CIM-FEM method to solve the 2D Case :
%                with zero Dirichlet bundary condition
%            
%Explainations:
%   this expriment results showing in  Fig.2
%***Author：
%           F.G. MA
%*****Finishi Data：
%           2026-10-22
%**************************************************************************
% PLOT_EXAMPLE3_CIM_FEM_CLG_ONECASE_T05
% Example 3: six-peak, two-dimensional homogeneous multi-term TFJE.
% Same CIM in time, P1 FEM vs Legendre--Galerkin with Chebyshev input (CLG).
% Plots ONE selected parameter case at t=0.25 (FEM and CLG solutions).
%
% Main outputs:
%   Example3_t05_onecase_CIM_FEM_contours.png   (one selected case)
%   Example3_t05_onecase_CIM_CLG_contours.png   (one selected case)
%   Example3_t05_onecase_FEM_vs_CLG_contours.png (FEM left, CLG right)
% Optional subsequent quantitative output:
%   Example3_t05_onecase_relative_L2_errors.png
%   Example3_t05_onecase_results.mat
%
% No additional MATLAB toolboxes are required (MATLAB R2025a compatible).
% MATLAB execution was not available to the code author; verify locally.

clc; close all;

%% ------------------ Settings ------------------------------------------
opt.M                = 14;       % exactly (M-1)^2 = 361 DOF per method
opt.N                = 50;       % same CIM nodes for FEM and CLG
opt.t                = 0.25;      % the prescribed observation time
opt.nplot            = 1000;      % common visualization grid
opt.doQuantitative   = false;    % only solution contours; set true for later L2 comparison
opt.Mref             = 128;       % high-order numerical reference, errors only
opt.Nref             = 100;
opt.showSeparate     = true;     % separate FEM and CLG single-panel figures
opt.saveFigures      = true;
opt.outputDir        = pwd;
opt.caseIndex        = 1;        % 1:(.25,.55,.75), 2:(.50,.35,.45), 3:(.75,.15,.15)

% Select exactly ONE of the three fractional-order cases in Example 3.
allOrderSets = [0.25 0.55 0.75;
                0.50 0.35 0.45;
                0.75 0.15 0.15];
assert(ismember(opt.caseIndex,1:3),'opt.caseIndex must be 1, 2 or 3.');
orderSets = allOrderSets(opt.caseIndex,:);
nCases = 1;
caseNames = {sprintf('(%c) $alpha$=%.2f, $beta$=%.2f, $gamma$=%.2f', ...
    'a'+opt.caseIndex-1, orderSets(1), orderSets(2), orderSets(3))};

par.a=1.0; par.b=0.5;
par.alpha1=0.10; par.alpha2=0.18;
par.beta1=0.05; par.beta2=0.10;
par.a1=0.35; par.a2=0.20;
par.b1=0.25; par.b2=0.15;

assert(opt.M>=4 && opt.N>=8 && opt.t>0);
assert(~opt.doQuantitative || (opt.Mref>opt.M && opt.Nref>=opt.N));
if opt.saveFigures && ~exist(opt.outputDir,'dir'); mkdir(opt.outputDir); end

nDof=(opt.M-1)^2;
fprintf('Example 3 | t=%.3f | N=%d | M=%d | %d spatial DOF each\n', ...
    opt.t,opt.N,opt.M,nDof);
fprintf('Both use the same CIM nodes; K=J=2, D=1, f=0.\n');

%% ------------------ Same CIM nodes and spatial DOF ----------------------
[z,dz,tau] = cim_hyperbola(opt.N,opt.t);
check = (tau/pi)*imag(sum(exp(z*opt.t).*dz./z));
assert(abs(check-1)<1e-7,'CIM 1/z inversion check failed.');
[MF,KF,bF,xF] = fem_p1_matrices(opt.M);
[MC,KC,bC]   = clg_matrices(opt.M);

xplot=linspace(-1,1,opt.nplot);
[XP,YP]=meshgrid(xplot,xplot); % rows correspond to y, columns to x
Phi=eval_dirichlet_legendre(xplot,opt.M);

UF=cell(1,nCases); UC=cell(1,nCases);
for j=1:nCases
    par.alpha=orderSets(j,1);
    par.beta=orderSets(j,2);
    par.gamma=orderSets(j,3);
    fprintf('Computing %s ...\n',caseNames{j});
    cF=cim_linear_solve(MF,KF,bF,z,dz,tau,opt.t,par);
    cC=cim_linear_solve(MC,KC,bC,z,dz,tau,opt.t,par);
    UF{j}=eval_fem_triangles(cF,xF,XP,YP);
    UC{j}=(Phi*reshape(cC,opt.M-1,opt.M-1)*Phi.').';
end

%% ------------------ First: plot the physical solutions -----------------
% IMPORTANT: no reference solution and no error metric in these contours.
% A shared colormap and shared limits make FEM and CLG directly comparable.
cMax=0;
for j=1:nCases
    cMax=max([cMax,max(UF{j}(:)),max(UC{j}(:))]);
end
cMax=max(cMax,eps);
% Values may vary with t: scale color limits to solution at t=0.5,
% not to the t=0.1 colorbar maximum in the original manuscript figure.
colorLimits=[0,cMax];
cmap=example3_colormap(256);

if opt.showSeparate
    hF=make_contour_figure('CIM--P1 FEM',UF,XP,YP,cmap,colorLimits, ...
        caseNames,opt,nDof);
    save_plot(hF,'Example3_t05_onecase_CIM_FEM_contours.png',opt);

    hC=make_contour_figure('CIM--CLG',UC,XP,YP,cmap,colorLimits, ...
        caseNames,opt,nDof);
    save_plot(hC,'Example3_t05_onecase_CIM_CLG_contours.png',opt);
end

hBoth=figure('Name','Example 3 at t=0.5: FEM and CLG solutions', ...
    'Color','w','Position',[45 35 1230 585]);
tlo=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    if row==1; data=UF; method='CIM--P1 FEM';
    else; data=UC; method='CIM--CLG'; end
    for j=1:nCases
        ax=nexttile(tlo);
        draw_one_contour(ax,XP,YP,data{j},cmap,colorLimits);
        title(ax,sprintf('%s | %s',method,caseNames{j}), ...
            'FontSize',10,'Interpreter','none');
    end
end
cb=colorbar;
cb.Layout.Tile='east';
cb.Label.String='p(t,x,y)';
%sgtitle(tlo,sprintf('Six-peak MT-TFJE at t=%.2f | %d DOF per method | N=%d', ...
%    opt.t,nDof,opt.N),'FontWeight','bold');
save_plot(hBoth,'Example3_t05_onecase_FEM_vs_CLG_contours.png',opt);
drawnow;

%% ------------------ Second: quantitative comparison --------------------
errF=[];errC=[];
if opt.doQuantitative
    fprintf('\nPhysical solution contours have been produced.\n');
    fprintf('Now computing HIGH-ORDER REFERENCE ONLY for the error comparison...\n');
    [MR,KR,bR]=clg_matrices(opt.Mref);
    [zr,dzr,taur]=cim_hyperbola(opt.Nref,opt.t);
    Pr=eval_dirichlet_legendre(xplot,opt.Mref);
    errF=zeros(nCases,1); errC=zeros(nCases,1);
    for j=1:nCases
        par.alpha=orderSets(j,1);
        par.beta=orderSets(j,2);
        par.gamma=orderSets(j,3);
        cR=cim_linear_solve(MR,KR,bR,zr,dzr,taur,opt.t,par);
        UR=(Pr*reshape(cR,opt.Mref-1,opt.Mref-1)*Pr.').';
        errF(j)=relative_l2(UF{j},UR,xplot);
        errC(j)=relative_l2(UC{j},UR,xplot);
        fprintf('Case %c at t=%.2f: FEM L2 = %.5e; CLG L2 = %.5e\n', ...
            'a'+opt.caseIndex-1,opt.t,errF(j),errC(j));
    end
    he=figure('Name','Relative L2 error (reference used only here)', ...
        'Color','w','Position',[150 130 680 470]);
    bar([errF,errC],'grouped');
    set(gca,'XTick',1,'XTickLabel',{sprintf('(%c)', 'a'+opt.caseIndex-1)});
    ylabel('Relative L^2 error against high-order CLG reference');
    xlabel('Fractional-order case');
    title(sprintf('Error comparison at t=%.2f, %d DOF, N=%d',opt.t,nDof,opt.N));
    legend('CIM--P1 FEM','CIM--CLG','Location','best');
    grid on;
    save_plot(he,'Example3_t05_onecase_relative_L2_errors.png',opt);
end

if opt.saveFigures
    save(fullfile(opt.outputDir,'Example3_t05_onecase_results.mat'), ...
        'opt','par','orderSets','xplot','UF','UC','errF','errC', '-v7');
end
fprintf('Completed. Exactly one parameter case was plotted.\n');
end

function fig=make_contour_figure(method,fields,XP,YP,cmap,limits,caseNames,opt,nDof)
fig=figure('Name',[method,' solution contours at t=0.5'], ...
    'Color','w','Position',[80 80 620 560]);
lo=tiledlayout(1,numel(fields),'TileSpacing','compact','Padding','compact');
for j=1:numel(fields)
    ax=nexttile(lo);
    draw_one_contour(ax,XP,YP,fields{j},cmap,limits);
    title(ax,caseNames{j},'FontSize',11,'Interpreter','none');
end
cb=colorbar;
cb.Layout.Tile='east';cb.Label.String='p(t,x,y)';
sgtitle(lo,sprintf('%s at t=%.2f | %d DOF | N=%d', ...
    method,opt.t,nDof,opt.N),'FontWeight','bold');
end

function draw_one_contour(ax,X,Y,U,cmap,limits)
% Match the paper's Fig. 5: a muted spectral/rainbow palette with fine contours.
% Not an imagesc heatmap and not an error or reference-solution plot.
Udraw=U;
if min(U(:)) < -1e-10
    warning('Negative solution values detected; check CIM accuracy.');
else
    Udraw=max(U,0); % suppress roundoff-only negatives in the plot
end
cMax=limits(2);
fillLevels=linspace(0,cMax,64);
lineLevels=linspace(0,cMax,33);
contourf(ax,X,Y,Udraw,fillLevels,'LineStyle','none');
hold(ax,'on');
contour(ax,X,Y,Udraw,lineLevels,'LineColor',[0.80 0.86 0.85], ...
    'LineWidth',0.24);
hold(ax,'off');
colormap(ax,cmap);
caxis(ax,limits);
axis(ax,'equal');
xlim(ax,[-1 1]);ylim(ax,[-1 1]);
set(ax,'YDir','normal');
xlabel(ax,'x');ylabel(ax,'y');
end

function cmap=example3_colormap(n)
% Approximate the ORIGINAL Example 3's muted blue/teal/green/yellow/red
% from the reference image supplied by the author (bottom -> top).
% Color endpoints are not taken from MATLAB's highly saturated turbo map.
loc=[0 .08 .17 .26 .35 .44 .53 .62 .70 .79 .88 1];
rgb=[53 60 94;
     43 67 125;
     48 97 157;
     51 130 177;
     58 158 181;
     75 183 179;
     118 202 174;
     176 217 165;
     235 227 156;
     249 221 133;
     250 175 99;
     196 80 71]/255;
cmap=interp1(loc,rgb,linspace(0,1,n),'pchip');
cmap=min(1,max(0,cmap));
end

function save_plot(fig,filename,opt)
if opt.saveFigures
    exportgraphics(fig,fullfile(opt.outputDir,filename),'Resolution',250);
end
end

%% ========================== Subfunctions ================================

function val=initial_six_peaks(x,y)
% EXACT initial condition from Example 3, Eq. (92), including pi/2.
s=zeros(size(x));
for k=0:5
    cx=0.8*cos(k*pi/3);
    cy=0.8*sin(k*pi/3);
    s=s+exp(-64*((x-cx).^2+(y-cy).^2));
end
val=(pi/2)*(1-x.^2).*(1-y.^2).*s;
end

function eta=laplace_symbol(z,p)
num=1+p.a*z^p.alpha+p.a1*z^p.alpha1+p.a2*z^p.alpha2;
den=1+p.b*z^p.beta+p.b1*z^p.beta1+p.b2*z^p.beta2;
eta=z^p.gamma*(num/den);
end

function [z,dz,tau]=cim_hyperbola(N,times)
% Paper's hyperbola (46) and optimized parameters (53),(55).
% alphaTilde=pi/4, deltaPrime=pi/8, 0<d<pi/2-alpha-deltaPrime.
at=pi/4; dp=pi/8;
d=0.9*min(at,pi/2-at-dp);
nu=0.5; machine_eps=1e-14;
t0=min(times); Lambda=max(times)/t0;
objective=@(rho) contour_objective(rho,N,Lambda,at,d,nu,machine_eps);
rho=fminbnd(objective,0.02,0.98); % MATLAB base functionality
ar=acosh(Lambda/((1-nu)*(1-rho)*sin(at-d)));
tau=ar/N;
mu=2*pi*d*N*(1-rho)/(t0*Lambda*ar);
phi=((0:N-1)+0.5)*tau;
z=mu*(1+sin(1i*phi-at));
dz=1i*mu*cos(1i*phi-at);
end

function v=contour_objective(rho,N,Lambda,at,d,nu,er)
ar=acosh(Lambda/((1-nu)*(1-rho)*sin(at-d)));
R=2*pi*d*N/ar;
v=er*exp((1-rho)*R)+exp(-rho*R)/(1-exp(-R));
end

function coeff=cim_linear_solve(Mass,Stiff,b0,z,dz,tau,times,par)
% One complex solve per contour node, shared across all requested times.
nt=numel(times); coeff=zeros(numel(b0),nt);
for k=1:numel(z)
    ek=laplace_symbol(z(k),par);
    A=Stiff+ek*Mass;
    rhs=(ek/z(k))*b0;    % transformed homogeneous initial-data load
    chat=A\rhs;          % sparse complex MATLAB direct solver
    timeWeight=(tau/pi)*exp(z(k)*times(:).')*dz(k);
    coeff=coeff+imag(chat*timeWeight);
end
end

function [Mass,Stiff,b0,x]=fem_p1_matrices(M)
% Uniform triangular mesh on [-1,1]^2, split each square A--C.
% m=M-1 interior points in each direction => m^2 unknowns.
x=linspace(-1,1,M+1);
[X,Y]=meshgrid(x,x);
vertices=[X(:),Y(:)];
ids=reshape(1:(M+1)^2,M+1,M+1);
triangles=zeros(2*M*M,3);
t=0;
for ix=1:M
    for iy=1:M
        A=ids(iy,ix); B=ids(iy,ix+1);
        C=ids(iy+1,ix+1); D=ids(iy+1,ix);
        triangles(t+1,:)=[A,B,C];
        triangles(t+2,:)=[A,C,D];
        t=t+2;
    end
end
nTri=size(triangles,1);nV=size(vertices,1);
I=zeros(9*nTri,1);J=I;SV=I;MV=I;
load0=zeros(nV,1);
% Seven-point Dunavant rule (weights sum to one, multiplied by area).
a=0.059715871789770; b=0.470142064105115;
c=0.797426985353087; d=0.101286507323456;
bary=[1/3,1/3,1/3;
      a,b,b; b,a,b; b,b,a;
      c,d,d; d,c,d; d,d,c];
qw=[0.225000000000000;
    0.132394152788506*ones(3,1);
    0.125939180544827*ones(3,1)];

for it=1:nTri
    idx=triangles(it,:).';P=vertices(idx,:);
    x1=P(1,1);x2=P(2,1);x3=P(3,1);
    y1=P(1,2);y2=P(2,2);y3=P(3,2);
    area=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2;
    grads=[y2-y3,y3-y1,y1-y2;
           x3-x2,x1-x3,x2-x1]/(2*area);
    ke=area*(grads.'*grads);
    me=(area/12)*[2 1 1;1 2 1;1 1 2];
    qpt=bary*P;
    fv=initial_six_peaks(qpt(:,1),qpt(:,2));
    le=area*(bary.'*(qw.*fv));
    load0(idx)=load0(idx)+le;
    pos=(it-1)*9+(1:9);
    I(pos)=repmat(idx(:),3,1);
    J(pos)=repelem(idx(:),3,1);
    SV(pos)=ke(:);
    MV(pos)=me(:);
end
Mass=sparse(I,J,MV,nV,nV);
Stiff=sparse(I,J,SV,nV,nV);
interior=ids(2:M,2:M);
interior=interior(:);
Mass=Mass(interior,interior);
Stiff=Stiff(interior,interior);
b0=load0(interior);
end

function [Mass,Stiff,b0]=clg_matrices(M)
% Dirichlet Legendre basis phi_j=c_j*(P_j-P_{j+2}), j=0,..,M-2.
% 1D stiffness is the identity; 1D mass is pentadiagonal.
m=M-1;
c=1./sqrt(4*(0:m-1)+6);
B=spalloc(m,m,3*m);
for j=0:m-1
    B(j+1,j+1)=c(j+1)^2*(2/(2*j+1)+2/(2*j+5));
    if j+2<=m-1
        v=-c(j+1)*c(j+3)*2/(2*j+5);
        B(j+1,j+3)=v;
        B(j+3,j+1)=v;
    end
end
E=speye(m);
Mass=kron(B,B);
Stiff=kron(B,E)+kron(E,B);

% CLG load: interpolate p0 on tensor-product CGL points, use FFT/DCT-I
% to recover Chebyshev coefficients, then project the interpolant exactly
% against the homogeneous Dirichlet Legendre test functions.
xc=cos(pi*(0:M)/M);
[Xc,Yc]=ndgrid(xc,xc);
fCGL=initial_six_peaks(Xc,Yc);
chebCoeffs=cheb2_coeff_fft(fCGL);

% Gauss-Legendre q=M+2 integrates products T_k*phi_j (degree <=2M)
% exactly to round-off. There is no spatial quadrature aliasing here.
[qx,qw]=gauss_legendre(M+2);
Phi=eval_dirichlet_legendre(qx,M);
T=cos(acos(qx(:))*(0:M));
Proj=Phi.'*(qw(:).*T); % m x (M+1)
B0=Proj*chebCoeffs*Proj.';
b0=B0(:);
end

function CC=cheb2_coeff_fft(F)
% FFT-based DCT-I on Lobatto values. Produces full Chebyshev coefficients
% F(x_j,y_l)=sum_{k,r=0}^M CC(k+1,r+1)*T_k(x_j)*T_r(y_l).
CC=cheb1_coeff_fft(F);
CC=cheb1_coeff_fft(CC.').';
end

function C=cheb1_coeff_fft(F)
M=size(F,1)-1;
Fext=[F; F(M:-1:2,:)];    % even extension to length 2M
D=real(fft(Fext,[],1))/M;
C=D(1:M+1,:);
C([1,M+1],:)=C([1,M+1],:)/2;
end

function [x,w]=gauss_legendre(n)
% Golub-Welsch: Gauss-Legendre quadrature on [-1,1].
k=1:n-1;
beta=k./sqrt(4*k.^2-1);
J=diag(beta,1)+diag(beta,-1);
[V,D]=eig(J);
[x,order]=sort(diag(D));
w=2*(V(1,order).').^2;
end

function Phi=eval_dirichlet_legendre(x,M)
% Size: numel(x) x (M-1)
x=x(:);
Leg=zeros(numel(x),M+1);
Leg(:,1)=1;
Leg(:,2)=x;
for k=2:M
    Leg(:,k+1)=((2*k-1)*x.*Leg(:,k)-(k-1)*Leg(:,k-1))/k;
end
m=M-1;
c=1./sqrt(4*(0:m-1)+6);
Phi=(Leg(:,1:m)-Leg(:,3:m+2)).*c;
end

function U=eval_fem_triangles(coeff,xnodes,X,Y)
% Evaluate P1 FEM solution exactly on the original A--C triangular grid.
M=numel(xnodes)-1;
nd=M+1;
Ugrid=zeros(nd,nd);
Ugrid(2:M,2:M)=reshape(coeff,M-1,M-1);
h=xnodes(2)-xnodes(1);
ix=min(M,max(1,floor((X-xnodes(1))/h)+1));
iy=min(M,max(1,floor((Y-xnodes(1))/h)+1));
xi=max(0,min(1,(X-xnodes(ix))/h));
eta=max(0,min(1,(Y-xnodes(iy))/h));
UA=Ugrid(sub2ind(size(Ugrid),iy,ix));
UB=Ugrid(sub2ind(size(Ugrid),iy,ix+1));
UC=Ugrid(sub2ind(size(Ugrid),iy+1,ix+1));
UD=Ugrid(sub2ind(size(Ugrid),iy+1,ix));
U=zeros(size(X));
lo=(eta<=xi);
U(lo)=(1-xi(lo)).*UA(lo)+(xi(lo)-eta(lo)).*UB(lo)+eta(lo).*UC(lo);
up=~lo;
U(up)=(1-eta(up)).*UA(up)+xi(up).*UC(up)+(eta(up)-xi(up)).*UD(up);
end

function err=relative_l2(U,Uref,x)
% Trapezoidal L2 norm on the SAME fine grid for both methods.
E=U-Uref;
num=trapz(x,trapz(x,abs(E).^2,2));
den=trapz(x,trapz(x,abs(Uref).^2,2));
err=sqrt(num/max(den,realmin));
end
