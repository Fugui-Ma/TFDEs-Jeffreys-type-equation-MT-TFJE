function [u]=Exact_solution(t,a,x)
u=t.^a*sin(pi*x);
end