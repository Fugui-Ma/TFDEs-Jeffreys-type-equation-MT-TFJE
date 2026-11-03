function err=J9_source_error(coeff,M,dim,caseName,p,t)
% J9_SOURCE_ERROR L2 error of the 1D/2D manufactured exact solutions.
% Performed ONCE after benchmarking, never inside the per-sample timer.
%
% Exact p(t,x)=t^kappa*v(x); p(t,x,y)=t^kappa*v(x)*v(y).
% Does not apply to the six-peak initial-value problem.

if strcmpi(caseName,'sixpeak')
    err=NaN;
    return;
end
m=M-1;
[x,w]=gaussleg(max(180,3*M+30));
P=legendreV(x,M);
j=0:M-2;
Phi=bsxfun(@rdivide,P(:,1:M-1)-P(:,3:M+1),sqrt(4*j+6));

if strcmpi(caseName,'analytic')
    v=(1-x.^2).*(1.1+sin(pi*x));
elseif strcmpi(caseName,'finite')
    v=(1-x.^2).*(.7+abs(x-.173).^2.6);
else
    error('Unknown exact solution type.');
end

if dim==1
    d=Phi*coeff(1,:).'-t^p.kappa*v;
    err=sqrt(sum(w.*abs(d).^2));
elseif dim==2
    actual=Phi*reshape(coeff(1,:),m,m)*Phi.';
    exact=t^p.kappa*(v*v.');
    d=actual-exact;
    err=sqrt(sum(sum((w*w.').*abs(d).^2)));
else
    error('dim must be 1 or 2');
end
end

function P=legendreV(x,M)
P=zeros(numel(x),M+1);
P(:,1)=1;
if M>=1,P(:,2)=x;end
for k=1:M-1
    P(:,k+2)=((2*k+1)*x.*P(:,k+1)-k*P(:,k))/(k+1);
end
end

function [x,w]=gaussleg(n)
j=(1:n-1)';
b=j./sqrt(4*j.^2-1);
J=diag(b,1)+diag(b,-1);
[V,D]=eig(J);
[x,id]=sort(diag(D));
V=V(:,id);
w=2*(V(1,:).').^2;
end
