function L2_error = CLG_L2_Error(u_coef,M,kappa,current_time)
%**************************************************************************
% L2 error evaluated on an independent, higher-order LGL grid.
% The numerical solution is reconstructed directly in the CLG basis.
%**************************************************************************
M_error = max(128,2*M+2);
[x_error,w_error] = LGL_points_Weights(M_error);

UU_error = zeros(size(x_error));
for k = 0:M-2
    phi_k = Legendre_polynomial(k,x_error) ...
          - Legendre_polynomial(k+2,x_error);
    UU_error = UU_error + u_coef(k+1)*phi_k;
end

Exact_solu_error = Exact_solution(current_time,kappa,x_error);
L2_error = sqrt(sum(w_error.*abs(UU_error-Exact_solu_error).^2));
end
