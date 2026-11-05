function [Mm,S] = CLG_Matrices(M)
%**************************************************************************
% Assemble the CLG mass and stiffness matrices.
% This is the same Legendre basis as in the uploaded original program:
% phi_k(x)=L_k(x)+a_k*L_{k+1}(x)+b_k*L_{k+2}(x),
% with a_k=0 and b_k=-1.
%**************************************************************************
ak = zeros(1,M-1);
bk = -ones(1,M-1);

S = zeros(M-1,M-1);
for k = 0:M-2
    for j = 0:M-2
        if k == j
            S(k+1,j+1) = -(4*k+6)*bk(k+1);
        end
    end
end

Mm = zeros(M-1,M-1);
for j = 0:M-2
    for k = 0:M-2
        if j == k
            Mm(j+1,k+1) = 2/(2*k+1) ...
                + ak(k+1)^2*2/(2*k+3) ...
                + bk(k+1)^2*2/(2*k+5);

        elseif j == k+1
            Mm(j+1,k+1) = ak(k+1)*2/(2*k+3) ...
                + ak(k+2)*2/(2*k+5);
            Mm(k+1,j+1) = Mm(j+1,k+1);

        elseif j == k+2
            Mm(j+1,k+1) = bk(k+1)*2/(2*k+5);
            Mm(k+1,j+1) = Mm(j+1,k+1);

        else
            Mm(j+1,k+1) = 0;
        end
    end
end
end
