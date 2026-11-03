%**************************************************************************
%**    Program£º Main code-By using the CGL-method to solve the 1D Case :
%                with zero Dirichlet bundary condition
%            
%Explainations:
%   this expriment results showing in  Table 1 and Table 3
%***Author£º
%           F.G. MA
%*****Finishi Data£º
%           2022-8-16
%**************************************************************************
close all, clear all, clc
% % parameters
% chi=0.17; 
% a=1/2;
% %Alpha=1; Beta=1; Gamma=1;
% Alpha=1; Beta=1; Gamma=1;
% parameters
kappa=1/4; a=1; b=1000; 

Alpha=0.25; Beta=0.55; Gamma=0.75;
%Alpha=0.50; Beta=0.35; Gamma=0.450;
%Alpha=0.75; Beta=0.15; Gamma=0.15;
%Alpha=1.00; Beta=1.00; Gamma=1.00 ;

Nmax=50; current_time=0.5;

M=[4,6,8,10,12,14,16];

for kk=1:length(M)
   [Err_inf(kk),Err_L2(kk)]=Space_error_order_CGL(kappa,a,b,Alpha,Beta,Gamma,Nmax,M(kk),current_time); 
end
for j=1:length(M)-1
       Ratio_L2(j)=log2(Err_L2(j)/Err_L2(j+1))/log2(M(j+1)/M(j));
       Ratio_Linf(j)=log2(Err_inf(j)/Err_inf(j+1))/log2(M(j+1)/M(j));
end

 format short e
disp('The L2 error for p in space are:')
   Err_L2
disp('The L_inf error for p in space are:')
   Err_inf
   
 format long
disp('The ratio of L2 error for p in space are:')
   Ratio_L2
disp('The ratio of L_inf error for p in space are:')
   Ratio_Linf

