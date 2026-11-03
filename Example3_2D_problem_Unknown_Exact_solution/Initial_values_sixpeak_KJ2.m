function f = Initial_values_sixpeak_KJ2(x,y)
% Vectorized six-Gaussian-peak initial condition on (-1,1)^2.
f = zeros(size(x+y));
for k = 0:5
    xk = 0.8*cos(k*pi/3);
    yk = 0.8*sin(k*pi/3);
    f = f + exp(-64*((x-xk).^2+(y-yk).^2));
end
f = pi.^2.*(1-x.^2).*(1-y.^2).*f;
end