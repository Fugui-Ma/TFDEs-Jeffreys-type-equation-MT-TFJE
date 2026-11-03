function [u]=Exact_solution(t,a,x)
%u=(t.^a)*(sin(pi*x)+1.1);
u=(t.^a)*(sin(pi*x)+1.1).*(1-x.^2);
end