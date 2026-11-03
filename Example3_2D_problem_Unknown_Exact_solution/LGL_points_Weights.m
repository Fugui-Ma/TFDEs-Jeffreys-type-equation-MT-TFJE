function [nodes,weights] = LGL_points_Weights(M)
% Legendre-Gauss-Lobatto points and quadrature weights on [-1,1].
% Interior nodes solve x*P_M(x)-P_{M-1}(x)=0 (P'_M(x)=0).
% Newton is initialized with Chebyshev-Lobatto points.
validateattributes(M,{'numeric'},{'scalar','integer','>=',4},mfilename,'M');
nodes = -cos(pi*(0:M)'/M);
interior = nodes(2:end-1);
for iter = 1:100
    [pM,pMm1] = legendre_pair(M,interior);
    correction = (interior.*pM-pMm1)./((M+1).*pM);
    interior = interior - correction;
    if max(abs(correction)) < 2e-14
        break;
    end
end
if iter == 100
    error('LGL node Newton iteration did not converge for M=%d.',M);
end
nodes = [-1;interior;1];
[nP,~] = legendre_pair(M,nodes);
weights = 2./(M*(M+1).*nP.^2);
end

function [pN,pNm1] = legendre_pair(N,x)
% Returns P_N(x) and P_{N-1}(x), using the stable three-term recursion.
pNm1 = ones(size(x));
pN = x;
for n = 2:N
    next = ((2*n-1).*x.*pN - (n-1).*pNm1)./n;
    pNm1 = pN;
    pN = next;
end
end
