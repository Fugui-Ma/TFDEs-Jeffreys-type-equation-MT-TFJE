function [R1,R2]=Main_Jeffreys_compare_all(showTransform)
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
% Run both examples. Each chart opens in a SEPARATE MATLAB figure window.
% Default: two windows (1D total runtime, 2D total runtime).
% Optional true: two additional transform-only windows.
if nargin<1,showTransform=false;end
R1=Main_Jeffreys_compare_1d(showTransform);
R2=Main_Jeffreys_compare_2d(showTransform);
end
