% PROGRAM : Main_code_computing the first internal states show the accuracy
%           of CIM 
%           This program is used to determin the numerical accuracy of the
%           CIM, where the optimal parameteres were choose in 
%  {LOPEZ FERNANDEZ, M., PALENCIA, C. & SCHADLE, A.(2006) A spectral order method for inverting sectorial
%  Laplace transforms. SIAM J. Numer. Anal., 44, 1332–1350 (electronic)， Ma.F.G., J. Sci., Comput.,45(2023)}
%DATA:  
%     ————2022-06-20————
%AUTHER：
%              F.G.Ma
%#########################################################################
clear all; clc

 N = 150;     % the number of the discrete points for the CIM 
% t0 = 0.1;   % start time 
 T=1;  
 t = 0.5; 
% B = 6;   % The time interval for the CIM % t=1.5
 alpha = 0.50;  beta = 0.35;  gama = 0.45;   % Fractional orders 
 
 %chi= 1.5;
 a=1; b=100;
 u_0 =1;       % initial value 

%  alpha1 = 0.6767; % Angle of the contour alpha
%  d = 0.7901;   % the width of the interval 
alpha1 = pi/4; % Angle of the contour alpha
delta_per=pi/8;
d = min(alpha1,pi/2-alpha1-delta_per); 

 
 spgridnum = [1,10,100,500];
  numlevel = length(spgridnum);             %统计字符长度.
  errorall = cell(1,numlevel);              %存放误差项
% time_CI = zeros(numlevel,1);              %存时间
 
  for j=1:numlevel
  % grid mesh;
    B =spgridnum(j)   % the CIM grid number is
    [Err,mu]=Hyperbola_contour_Jeff(N,t,T,alpha,beta,gama,a,b,d,alpha1,B,u_0);                       
    errorall(1,j) = {Err};    
  end
  
figure(1)
figure('Position', [100, 100, 500, 500]); % 设置图形大小
   semilogy(1:4:N,errorall{1,1}(1:4:end),'->',1:4:N,errorall{1,2}(1:4:end),'-o',...
            1:4:N,errorall{1,3}(1:4:end),'-d',1:4:N,errorall{1,4}(1:4:end),'-s','lineWidth',1.5); 
   h1 = xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 18);
   h2 = ylabel('$\mathrm{Error}(N)$ (log scale)', 'Interpreter', 'latex', 'FontSize', 18);
   legend('\Lambda = 1','\Lambda = 10','\Lambda = 100','\Lambda = 500');
  
   hold off
   
%   hold on
   
%    for k=[1,5,10,50]
%    plot(20:5:60,exp(20-4.4./(log(20:5:60)+log(k)).*(20:5:60)),'-')
%    hold on
%    end
%    hold off

 %#######################################################################
   
 figure(2)  
 figure('Position', [100, 100, 500, 500]); % 设置图形大小
   paras=[0.25;0.5;0.75;1]';
   numlevel = length(paras);
   for j=1:numlevel 
       paras1=paras(j);
       [Err,mu]=Hyperbola_contour_Jeff(N,t,T,paras1,paras1,paras1,a,b,d,alpha1,10,u_0);                         
      erroralll(1,j) = {Err};   
   end 
   paras3=[0.2;0.4;0.6;0.8]';
   numlevel = length(paras3);
   for j=1:numlevel 
       paras2=paras3(j);
       [Err,mu]=Hyperbola_contour_Jeff(N,t,T,0.95-paras2,paras2,paras2+0.1*(j),a,b,d,alpha1,10,u_0);                         
      errorall3(1,j) = {Err};
   end 
 
  
   semilogy(1:4:N,erroralll{1,1}(1:4:end),'-d',1:4:N,erroralll{1,2}(1:4:end),'>',...
            1:4:N,erroralll{1,3}(1:4:end),'*',...
            1:4:N,errorall3{1,1}(1:4:end),'<',1:4:N,errorall3{1,2}(1:4:end),'x',...
            1:4:N,errorall3{1,3}(1:4:end),'^',1:4:N,errorall3{1,4}(1:4:end),'s','lineWidth',1.5); 
   xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 18);
   ylabel('$\mathrm{Error}(N)$ (log scale)', 'Interpreter', 'latex', 'FontSize', 18);
   hold off
   legend('\alpha = \beta = \gamma = 0.25','\alpha = \beta = \gamma = 0.50',...
          '\alpha = \beta = \gamma = 0.75',...
          '\alpha =0.75, \beta =0.20, \gamma = 0.21','\alpha =0.55, \beta =0.40, \gamma = 0.42',...
          '\alpha =0.35, \beta =0.60, \gamma = 0.63','\alpha =0.15, \beta =0.80, \gamma = 0.84')
  
   
 figure(3)  
 figure('Position', [100, 100, 500, 500]); % 设置图形大小

   a=1;
   b=[1,10,100,1000]
   numlevel = length(b);
   for j=1:numlevel 
       b1=b(j);
       [Err,mu]=Hyperbola_contour_Jeff(N,t,T,alpha,beta,gama,a,b1,d,alpha1,50,u_0);                         
      erroralll(1,j) = {Err};   
   end 
   semilogy(1:4:N,erroralll{1,1}(1:4:end),'-d',1:4:N,erroralll{1,2}(1:4:end),'>',...
            1:4:N,erroralll{1,3}(1:4:end),'*',1:4:N,erroralll{1,4}(1:4:end),'o','lineWidth',1.5);        
   xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 18);
   ylabel('$\mathrm{Error}(N)$ (log scale)', 'Interpreter', 'latex', 'FontSize', 18);
   hold off
   legend('$b/a=1$','$b/a=10$','$b/a=100$','$b/a=1000$')
          
   
   
   
   
   
   
   
   
   
   
    
   

 