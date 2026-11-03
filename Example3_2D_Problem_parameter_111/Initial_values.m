function f=Initial_values(x,y)
% Program: Initial values for 1D
%  Author: F.G. MA
%

%f=(0.5+0.3.*sin(pi.*x)).^(-1).*sin(pi.*y).*exp(-2.5*x.*y);

%f= ((0.5+0.3*cos(pi.*x)).^(-1).*(0.5+0.3*cos(pi.*y).^(-1))).*exp(-10*x.*y);

%f= (1./pi).^3.*((0.75+0.3*cos(pi.*x)).^(-1).*(0.75+0.3*cos(pi.*y).^(-1))).*exp(-10*x.*y);

%f=1/(4*pi)*((0.15 + 0.3./sin(pi.*x)).*(0.15 + 0.3./sin(pi.*y)));

%f = -pi^(-3) .* ...
%    (1-x.^2).^3 .* (1-y.^2).^3 .* ...
%    (0.75+0.3*cos(pi*x)).^(-1) .* ...
%    (0.75+0.3*sin(pi*y)).^(-1) .* ...
%    exp(-10*x.*y);

%f = cos(pi*x./2).*cos(pi*y./2).*(1-4/5*cos(pi.*x).*cos(pi.*y));

f = zeros(size(x+y));
for k = 0:5
    xk = 0.8*cos(k*pi/3);
    yk = 0.8*sin(k*pi/3);
    f = f + exp(-64*((x-xk).^2+(y-yk).^2));
end
f = pi.^2.*(1-x.^2).*(1-y.^2).*f;
end