function [f]=source_term(kappa,a,b,Alpha,Beta,Gamma,x,z)
% Program: source term for 1D
%  Author: F.G. MA
%
m=length(x);
n=length(z);
f=zeros(m,n);
for k=1:m
    for j=1:n
     %f(k,j)=(z(j).^(Alpha+Gamma-a-1)+z(j).^(Alpha-a-1)+pi.^2*z(j).^(-a-1)+...
     %       pi.^2*chi.^Beta*z(j).^(Beta-a-1))*gamma(a+1)*sin(pi*x(k));
     %f(k,j)=gamma(kappa+1)*((sin(pi*x(k))+1.1)*(z(j).^(Gamma-kappa-1)+a*z(j).^(Alpha+Gamma-kappa-1))+...
     %       sin(pi*x(k))*pi.^2*(z(j).^(-kappa-1)+b*z(j).^(Beta-kappa-1)));
     u = (sin(pi*x(k))+1.1).*(1-x(k).^2);
     uu = -sin(pi*x(k))*(pi.^2*(1-x(k).^2)+2)-2.2-4*pi*x(k)*cos(pi*x(k));
     f(k,j)=gamma(kappa+1)*(u*(z(j).^(Gamma-kappa-1)+a*z(j).^(Alpha+Gamma-kappa-1))-...
            uu*(z(j).^(-kappa-1)+b*z(j).^(Beta-kappa-1)));
    end
end
% f=gamma(a+1)*((z.^(Alpha+Gamma-a-1)+z.^(Alpha-a-1))...
%     *sin(pi.*x)-pi.^2*(z.^(-a-1)+chi.^Beta*z^(Beta-a-1))*sin(pi*x));
end
