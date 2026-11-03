function err = Jeffreys2D_DoublingError_KJ2(low,hi,spHi,xDense)
% Compare p_M^N and p_{2M}^N for the SAME N and t.
% The Shen basis phi_k=P_k-P_{k+2} is independent of M; therefore the
% coefficients of the low-order approximation can be padded with zeros.
% The continuous 2D L2 norm is exact via the tensor Galerkin mass matrix:
% ||D||^2 = trace(D'*Mm*D*Mm), up to floating-point roundoff.
% L-infinity is the sampled max on an explicitly recorded common grid.
if ~isfield(low,'Ucoef') || ~isfield(hi,'Ucoef')
    error('Need full Galerkin coefficients Ucoef in both solutions.');
end
mLo=size(low.Ucoef,1);
mHi=size(hi.Ucoef,1);
if mHi<=mLo || size(low.Ucoef,2)~=mLo || size(hi.Ucoef,2)~=mHi
    error('Invalid coefficient shapes or high/low polynomial degrees.');
end
if spHi.M ~= mHi+1
    error('Mass matrix does not match higher polynomial degree.');
end
D=hi.Ucoef;
D(1:mLo,1:mLo)=D(1:mLo,1:mLo)-low.Ucoef;
Mm=full(spHi.Mm);
err.L2=sqrt(max(0,real(sum(sum(conj(D).*(Mm*D*Mm))))));
xDense=xDense(:);
P=Legendre_Vandermonde_Jeffreys(xDense,spHi.M);
B=P(:,1:mHi)-P(:,3:mHi+2);
E=B*D*B.';
err.Linf_sample=max(abs(E(:)));
err.gridCount=numel(xDense);
err.normDescription='Continuous L2 exact via Galerkin mass; Linf sampled grid maximum.';
end
