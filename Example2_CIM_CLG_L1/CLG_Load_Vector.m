function F = CLG_Load_Vector(f,xx,M)
%**************************************************************************
% Backward Legendre transform and CLG projection of a nodal source vector.
% This follows the transform in the uploaded original CIM program.
%**************************************************************************
g_hat = complex(zeros(M+1,1));

for k = 0:M-1
    temp = 0;
    for j = 0:M
        L_k = Legendre_polynomial(k,xx(j+1));
        L_M = Legendre_polynomial(M,xx(j+1));
        temp = temp + f(j+1)*L_k/(L_M)^2;
    end
    g_hat(k+1) = (2*k+1)/(M*(M+1))*temp;
end

temp = 0;
for j = 0:M
    L_M = Legendre_polynomial(M,xx(j+1));
    temp = temp + f(j+1)/L_M;
end
g_hat(M+1) = 1/(M+1)*temp;

% Since a_k=0 and b_k=-1 in the original CLG basis:
F = zeros(M-1,1);
for k = 0:M-2
    F(k+1) = g_hat(k+1)/(k+1/2) ...
           - g_hat(k+3)/(k+5/2);
end
end
