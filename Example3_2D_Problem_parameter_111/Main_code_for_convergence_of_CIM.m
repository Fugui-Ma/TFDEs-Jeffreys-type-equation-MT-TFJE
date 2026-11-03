%**************************************************************************
%**    Program£º Main code-By using the CGL-method to solve the 2-D Case :
%                with zero Dirichlet bundary condition for CIM
%            
%Explainations:
%   this expriment results showing in  Table 4 
%***Author£º
%           F.G. MA
%*****Finishi Data£º
%           2022-11-10
%**************************************************************************
close all, clear all, clc
% chi=1;
% Alpha=1; Beta=1; Gamma=1;
 M=12; Nmax=50; current_time=0.5;

a=1; b=0.5; 

%Alpha=0.25; Beta=0.55; Gamma=0.75;
%Alpha=0.50; Beta=0.35; Gamma=0.45;
%Alpha=0.75; Beta=0.15; Gamma=0.15;
Alpha=1.00; Beta=1.00; Gamma=1.00 ;


UU = Code_2D_solver_zero_Initial_Smooth_CIM(a,b,Alpha,Beta,Gamma,M,Nmax,current_time);

UU_ref = Code_2D_solver_zero_Initial_Smooth_CIM(a,b,Alpha,Beta,Gamma,M,100,current_time);

for k=1:Nmax
    error_N(k)=norm(UU{k}(2:end-1)-UU_ref{end}(2:end-1));
end
format short e
[error_N(4),error_N(6),error_N(10),error_N(20),error_N(30),error_N(50)]
