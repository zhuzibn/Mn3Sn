function check_anisotropy_stability(source,reference,out)
% Analytic constrained Hessian of E/mu for the half-strength triangle.
% Every candidate patch is a sum of these energies; graph coverage is checked
% independently by anisotropy_trial.py. This establishes local, not global, stability.
constantfile;data=load(fullfile(source,'2x2-static.mat'));
A=.5*(data.A~=0);C=.5*sign(data.C);axes=data.axes;
assert(isequal(A,A') && isequal(C,-C') && size(A,1)==3);
mu=3*mub;J=100.812;D=.833e-3*ele/mu;K=2*.196e-3*ele/mu;
ki=K*sum(A~=0,2)/4;Q=zeros(9);X=[0 0 0;0 0 1;0 -1 0];
for i=1:3
 ii=(i-1)*3+(1:3);
 Q(ii,ii)=-ki(i)*(axes(i,:)'*axes(i,:));
 for j=1:3
  if i==j,continue;end
  jj=(j-1)*3+(1:3);Q(ii,jj)=J*A(i,j)*eye(3)-D*C(i,j)*X;
 end
end
assert(max(abs(Q-Q'),[],'all')<1e-12);
phis=30:60:330;spectra=zeros(6);residual=zeros(1,6);
for p=1:6
 f=fopen(fullfile(reference,sprintf('2x2-current-phi%d-p0-final.bin',phis(p))),'r','ieee-le');assert(f>=0);
 flat=fread(f,Inf,'double');fclose(f);assert(numel(flat)==9);m=reshape(flat,3,3)';
 B=reshape(-Q*flat,3,3)';P=zeros(9,6);L=Q;
 for i=1:3
  ii=(i-1)*3+(1:3);jj=(i-1)*2+(1:2);
  u=cross(m(i,:),[0 1 0]);u=u/norm(u);v=cross(m(i,:),u);
  P(ii,jj)=[u' v'];L(ii,ii)=L(ii,ii)+dot(m(i,:),B(i,:))*eye(3);
 end
 H=P'*L*P;assert(max(abs(H-H'),[],'all')<1e-12);
 spectra(p,:)=sort(eig((H+H')/2))';residual(p)=max(vecnorm(cross(m,B,2),2,2));
 assert(min(spectra(p,:))>1e-5 && residual(p)<1e-9);
end
result=struct('phi_deg',phis,'eigenvalues_T',spectra,'minimum_eigenvalue_meV',min(spectra,[],'all')*mu/(ele*1e-3),...
 'reference_torque_T',residual,'convention','Constrained second derivative of E/mu with orthonormal spin tangents');
f=fopen(fullfile(out,'stability.json'),'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(result,PrettyPrint=true));fclose(f);
fprintf('All six constrained triangle Hessians positive: min %.12g T (%.12g meV).\n',min(spectra,[],'all'),result.minimum_eigenvalue_meV);
end
