function prepare_driven(source,out)
% Actual-source fields, Gilbert identity and pulse-edge GPU fixtures.
% Run driven_figures.py prepare first; never calls main.m or saves root results.
root=pwd;addpath(root);constantfile;dev=gpuDevice(1);src=fileread('integrate_llg.m');
assert(~isfile(fullfile(out,'gpu-gates.json')),'Use a fresh output directory.');
writejson(fullfile(out,'environment.json'),struct('matlab',version,'gpu',dev.Name,'toolboxes',{ver}));
tags={'2x2','5x5','13x14','21x22'};checks={};
for q=1:numel(tags)
 tag=tags{q};dims=sscanf(tag,'%dx%d');natomW=dims(1);natomL=dims(2);natom=natomW*natomL;
 dwcalc=0;load_fixed_atom_distrib=0;save_fixed_atom_distrib=0;initial_theta_deg=210;
 systemgeneration;occ=find(atomtype_==1);N=numel(occ);bc=1;
 data=load(fullfile(source,[tag '-static.mat']));
 axes=[coeffs_1(occ) coeffs_2(occ) coeffs_3(occ)];assert(isequal(occ,data.occ) && isequal(axes,data.axes));
 modes={'triangle'};if N==3,modes{end+1}='current';end
 for modeindex=1:numel(modes)
  mode=modes{modeindex};A=.5*(data.A~=0);C=.5*sign(data.C);degree=sum(A~=0,2);
  assert(isequal(A,A') && isequal(C,-C') && all(degree==2 | degree==4));
  Jintra=100.812*(3*mub);Dsim=.833e-3*ele;musRE=3*mub;gamRE=gam;alp=.003;
  Ksim=zeros(natomW,natomL);Ksim(occ)=.196e-3*ele*degree/4;
  hex_constant(:)=.5;hdmi_constant(:)=.5;
  if strcmp(mode,'current')
   A=data.A;C=data.C;Ksim(occ)=.196e-3*ele;hex_constant(:)=1;hdmi_constant(:)=1;
  end
  rng(20261103);r=randn(N,3);r=r./vecnorm(r,2,2);
  names={'nextL','previousL','nextW','previousW','previousWleft','previousWright','nextWleft','nextWright'};
  for k=1:8,eval(['Asim' names{k} '=zeros(natomW,natomL);']);end
  eval(extract(src,'AsimnextL(2:2:natomW,:)=Jintra;','for ctL=1:natomL'));
  muigpu=ones(natomW,natomL);muigpu(occ)=musRE;
  mmxtmp=zeros(natomW,natomL);mmytmp=mmxtmp;mmztmp=mmxtmp;
  mmxtmp(occ)=r(:,1);mmytmp(occ)=r(:,2);mmztmp(occ)=r(:,3);
  eval(extract(src,'        mmxnextL=circshift','        if thermalenable'));
  actual=[hex_x(occ)+hani_x(occ),hex_y(occ)+hani_y(occ)+hdmi_y(occ),hex_z(occ)+hani_z(occ)+hdmi_z(occ)];
  predicted=-Jintra/musRE*A*r+Dsim/musRE*[zeros(N,1) C*r(:,3) -C*r(:,2)]+2*Ksim(occ)/musRE.*sum(r.*axes,2).*axes;
  ferr=max(abs(actual-predicted),[],'all');assert(ferr<1e-11);
  for trial=1:4
   signs=[1 1;1 -1;-1 -1;-1 1];jc=signs(trial,1)*5e10;
   Hext=[0 signs(trial,2)*.01 0];if trial==4,Hext(2)=.1;end
   BDSOTRE=hbar/2*.06*jc/(441850*30e-9);BDSTTRE=0;chi=0;
   field=actual+Hext;derivative=zeros(N,3);polar=repmat([1 0 0],N,1);
   for j=1:N
    [derivative(j,1),derivative(j,2),derivative(j,3)]=atom_llgs_rhs(r(j,1),r(j,2),r(j,3),1,0,0,0,0,1,gam/(1+alp^2),alp,field(j,1),field(j,2),field(j,3),BDSOTRE,0,0,0);
   end
   implicit=-gam*cross(r,field,2)+alp*cross(r,derivative,2)-gam*BDSOTRE*cross(r,cross(r,polar,2),2);
   residual=max(abs(implicit-derivative),[],'all')/gam;assert(residual<1e-10);
   m=data.seed+.01*r;m=m./vecnorm(m,2,2);
   mx_init(occ)=m(:,1);my_init(occ)=m(:,2);mz_init(occ)=m(:,3);
   dipolemode=0;thermalenable=0;enablefixedge=0;rk4=1;
   psjSHEx=1;psjSHEy=0;psjSHEz=0;psjSTTx=0;psjSTTy=0;psjSTTz=1;
   tstep=5e-15;gpusteps=11;gpusave=10*tstep;runtime=20*tstep;gpurun_number=2;savetstep=1;
   pulse_start=5*tstep;pulse_end=15*tstep;eval(src);
   mask=repmat(atomtype_==1,1,1,21);empty=~mask;norms=sqrt(mmx.^2+mmy.^2+mmz.^2);
   assert(max(abs(norms(mask)-1))<1e-12 && all(isfinite(norms),'all'));
   assert(all(mmx(empty)==0 & mmy(empty)==0 & mmz(empty)==0));
   assert(actual_updates==20 && pulse_active_updates==10 && numel(t)==21 && max(abs(t-(0:20)*tstep))<1e-25);
   trajectory=zeros(21,3*N);
   for j=1:21
    x=mmx(:,:,j);y=mmy(:,:,j);zspin=mmz(:,:,j);trajectory(j,:)=reshape([x(occ) y(occ) zspin(occ)]',1,[]);
   end
   name=sprintf('%s-%s-gpu%d',tag,mode,trial);
   f=fopen(fullfile(out,[name '.bin']),'w','ieee-le');assert(f>=0);fwrite(f,trajectory','double');fclose(f);
   checks{end+1}=struct('name',name,'size',tag,'mode',mode,'N',N,'field_y_T',Hext(2),'sot_T',BDSOTRE,'field_error_T',ferr,'gilbert_residual_over_gamma_T',residual,'updates',20,'active_updates',10,'vacancies_exact_zero',true); %#ok<AGROW>
  end
  fprintf('Driven source/GPU fixtures passed: %s %s, field error %.3g T\n',tag,mode,ferr);
 end
end
writejson(fullfile(out,'gpu-gates.json'),[checks{:}]);
end
function v=extract(s,a,b)
i=strfind(s,a);j=strfind(s,b);assert(numel(i)==1);j=j(find(j>i,1));v=s(i:j-1);
end
function writejson(path,data)
f=fopen(path,'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(data,PrettyPrint=true));fclose(f);
end
