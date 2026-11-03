function Verify_Jeffreys2D_KJ2_2()
% Validate (i) real function L2/Linf, (ii) batched modal CIM against the
% ORIGINAL tensor solver, and (iii) both historical order definitions.
sp=Build_Jeffreys2D_Spatial(12,@Initial_values_sixpeak_KJ2);
coefs=zeros(11,11);
coefs(1,1)=1;   % u(x,y)=(1-P2(x))(1-P2(y))
hi=struct('Ucoef',coefs);
lo=struct('Ucoef',zeros(5,5));
grid=linspace(-1,1,201).';
e=Jeffreys2D_DoublingError_KJ2(lo,hi,sp,grid);
assert(abs(e.L2-12/5)<1e-12,'Test L2=12/5 failed.');
assert(abs(e.Linf_sample-9/4)<1e-12,'Test Linf=9/4 failed.');
fprintf('PASS exact 2D L2, sampled Linf: %.12g, %.12g\n',e.L2,e.Linf_sample);

q=struct('K',2,'J',2,'a',1,'b',0.5, ...
  'Alpha',0.25,'Beta',0.55,'Gamma',0.75, ...
  'alpha_terms',[0.10 0.18],'beta_terms',[0.05 0.10], ...
  'a_terms',[0.35 0.20],'b_terms',[0.25 0.15], ...
  'N',16,'t',0.5,'t0',0.01,'B_gamma',150, ...
  'caseName','a','initial_function','Initial_values_sixpeak_KJ2');
times=[0.5,1.0];
modal=Jeffreys2D_CIM_ModalTimes_KJ2(sp,q,q.N,times);
Cmodal=reshape(modal(:,1),11,11);
C=sp.V*Cmodal*sp.V.';
r=Solve_Jeffreys2D_KJ2(sp,q);
assert(max(abs(C(:)-r.Ucoef(:)))<1e-11,'Time-batch CIM differs from Solve.');

q.N=10;
r10=Solve_Jeffreys2D_KJ2(sp,q);
m10=Jeffreys2D_CIM_ModalTimes_KJ2(sp,q,10,0.5);
q.N=100;
r100=Solve_Jeffreys2D_KJ2(sp,q);
m100=Jeffreys2D_CIM_ModalTimes_KJ2(sp,q,100,0.5);
modalError=norm(m100-m10,2);
D=r100.Ucoef-r10.Ucoef;
massError=sqrt(sum(sum(D.*(full(sp.Mm)*D*full(sp.Mm)))));
assert(abs(modalError-massError)<1e-11,'L2 modal orthonormality failed.');
assert(abs(log(16)/log(4)-2)<1e-14);
assert(abs(log(16)/log(2)-4)<1e-14);
fprintf('PASS same-Galerkin-space temporal error: %.8e (modal/mass).\n',modalError);
fprintf('PASS historical log4 and standard log2 order calculations.\n');
end
