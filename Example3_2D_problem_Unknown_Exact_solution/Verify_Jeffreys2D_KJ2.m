function Verify_Jeffreys2D_KJ2()
% Verifies tensor modal solver against the original Kronecker system for
% all three parameter groups, and the K=J=0 limit of the Laplace symbol.
M=6;
sp=Build_Jeffreys2D_Spatial(M,@Initial_values_sixpeak_KJ2);
p=struct('K',2,'J',2,'a',1,'b',.5, ...
     'Alpha',.25,'Beta',.55,'Gamma',.75, ...
     'alpha_terms',[.10,.18],'beta_terms',[.05,.10], ...
     'a_terms',[.35,.20],'b_terms',[.25,.15]);
triples=[.25,.55,.75; .50,.35,.45; .75,.15,.15];
z=1.3+0.7i;
Mm=sp.Mm;S=sp.S;
Mk=kron(Mm,Mm);
A=kron(S,Mm)+kron(Mm,S);
lam_sum=bsxfun(@plus,sp.lambda,sp.lambda.');
for ic=1:3
    p.Alpha=triples(ic,1);p.Beta=triples(ic,2);p.Gamma=triples(ic,3);
    [P,Q]=Jeffreys2D_PQ(z,p);
    lhs=z^p.Gamma*P*Mk+Q*A;
    rhs=z^(p.Gamma-1)*P*sp.F(:);
    u_direct=lhs\rhs;
    modal=sp.Fmode*(z^(p.Gamma-1)*P)./(z^p.Gamma*P+Q*lam_sum);
    Cfast=sp.V*modal*sp.V.';
    rel_error=norm(u_direct-Cfast(:),inf)/max(1,norm(u_direct,inf));
    assert(rel_error<1e-9, ...
         'Case (%d): tensor and original Kronecker systems disagree.',ic);
    fprintf('PASS: case (%s) direct vs tensor rel. error %.3e.\n', ...
          char('a'+ic-1),rel_error);
end
p0=p;
p0.a_terms=[0,0];p0.b_terms=[0,0];
[P0,Q0]=Jeffreys2D_PQ(z,p0);
errP=abs(P0-(1+p.a*z^p.Alpha));
errQ=abs(Q0-(1+p.b*z^p.Beta));
assert(errP<1e-13 && errQ<1e-13,'K=J=0 limit failed.');
fprintf('PASS: K=J=0 limiting symbol.\n');
end
