function [algebraicRate, exponentialDecay] = Jeffreys2D_ObservedOrders_KJ2(levels,E)
% Jeffreys2D_ObservedOrders_KJ2 - safe observational convergence diagnostics.
%
% levels: increasing positive row or column vector of M or N.
% E: nLevels x nCases array of positive errors.
% algebraicRate(i,:) = log(E(i-1,:)/E(i,:)) / log(levels(i)/levels(i-1))
% exponentialDecay(i,:) = log(E(i-1,:)/E(i,:)) / (levels(i)-levels(i-1))
% The first row is NaN. Near-roundoff errors produce NaN rather than
% meaningless giant rates. No manipulation of numerical solution values.
levels = levels(:);
if any(~isfinite(levels)) || any(levels<=0) || any(diff(levels)<=0)
    error('Levels must be finite, positive and strictly increasing.');
end
if size(E,1)~=numel(levels)
    error('E must have nLevels rows and one column per parameter case.');
end
n=numel(levels);
algebraicRate=nan(size(E));
exponentialDecay=nan(size(E));
for j=1:size(E,2)
    col=E(:,j);
    valid=col(isfinite(col) & col>0);
    if isempty(valid), continue; end
    noiseFloor=max(1e-12, 500*eps(max(1,max(valid))));
    for i=2:n
        eprev=col(i-1);
        enow=col(i);
        if isfinite(eprev) && isfinite(enow) && ...
                eprev>noiseFloor && enow>noiseFloor
            ratio=log(eprev/enow);
            algebraicRate(i,j)=ratio/log(levels(i)/levels(i-1));
            exponentialDecay(i,j)=ratio/(levels(i)-levels(i-1));
        end
    end
end
end
