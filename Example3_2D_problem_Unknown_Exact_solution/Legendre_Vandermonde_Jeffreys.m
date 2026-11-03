function L = Legendre_Vandermonde_Jeffreys(nodes,M)
% L(i,j+1) = P_j(nodes(i)), 0<=j<=M.
% Avoids dependence on Symbolic Math Toolbox / legendreP.
nodes = nodes(:);
n = numel(nodes);
L = zeros(n,M+1);
L(:,1) = 1;
if M>=1
    L(:,2) = nodes;
end
for k=2:M
    L(:,k+1) = ((2*k-1).*nodes.*L(:,k)-(k-1).*L(:,k-1))/k;
end
end
