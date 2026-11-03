function w=LGL_weights(n,x)

% w=legslbwt(n,x) returns n Legendre-Gauss-Lobatto weights, where the argument 
%   x is n LGL points (can be computed by x=legslb(n) or x=legslbndm(n)
% Here, we used the algorithm by ``E. Yakimiw, Accurate Computation of
% Weights in Classical Gauss-Christoffel Quadrature Rules, JCP, 129,
% 406-430, 1996. 
%  A more stable way to compute the quadrature weights
% Last modified on August 30, 2011

nb=n*(n-1); z=x(2:end-1);
[dy,y]=lepoly(n-1,z); 
f0=dy./y;
d0=(1-z.^2).^4; 
d2=nb*(1-z.^2).^3.*f0.^2/2; 
d3=-2*(1-z.^2).^2*nb.*z.*f0.^3/6; 
d4=nb*(1+3*z.^2+1.5*nb*(1-z.^2)).*(1-z.^2).*f0.^4/12;
d5=-nb*z.*(1+z.^2-7/6*nb*(1-z.^2)).*f0.^5/5;

w=2*(1-z.^2).^8./(nb*y.^2.*(d0+d2+d3+d4+d5).^2);

w=[2/nb;w;2/nb];

return;
  

