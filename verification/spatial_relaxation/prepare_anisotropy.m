function prepare_anisotropy(source,out)
% Isolated candidate: unchanged source loop, reciprocal bonds and K_i=K*z_i/4.
% prepare_spatial must previously have exported the source-derived matrices.
root=pwd;addpath(root);constantfile;dev=gpuDevice(1);src=fileread('integrate_llg.m');
meta=struct('matlab',version,'gpu',dev.Name,'toolboxes',{ver});
writejson(fullfile(out,'environment.json'),meta);
files=dir(fullfile(source,'*.model'));assert(numel(files)==6);
for q=1:numel(files)
 tag=erase(files(q).name,'.model');dims=sscanf(tag,'%dx%d');
 natomW=dims(1);natomL=dims(2);natom=natomW*natomL;
 dwcalc=0;load_fixed_atom_distrib=0;save_fixed_atom_distrib=0;initial_theta_deg=121;
 bc=1;systemgeneration;occ=find(atomtype_==1);N=numel(occ);
 data=load(fullfile(source,[tag '-static.mat']));
 axes=[coeffs_1(occ) coeffs_2(occ) coeffs_3(occ)];
 assert(isequal(occ,data.occ) && isequal(axes,data.axes));
 A=.5*(data.A~=0);C=.5*sign(data.C);z=sum(A~=0,2);
 assert(isequal(A,A') && isequal(C,-C') && all(z==2 | z==4));
 Jintra=100.812*(3*mub);Dsim=.833e-3*ele;musRE=3*mub;gamRE=gam;alp=.003;
 Ksim=zeros(natomW,natomL);Ksim(occ)=.196e-3*ele*z/4;
 hex_constant(:)=.5;hdmi_constant(:)=.5;
 rng(20261102);r=randn(N,3);r=r./vecnorm(r,2,2);
 % Execute the actual static source expressions with candidate coefficients.
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
 for trial=0:1
  m=data.seed;if trial,m=m+.01*r;m=m./vecnorm(m,2,2);end
  mx_init(occ)=m(:,1);my_init(occ)=m(:,2);mz_init(occ)=m(:,3);
  BDSOTRE=0;BDSTTRE=0;chi=0;dipolemode=0;thermalenable=0;enablefixedge=0;rk4=1;
  psjSHEx=1;psjSHEy=0;psjSHEz=0;psjSTTx=0;psjSTTy=0;psjSTTz=1;
  Hext=[0 0 0];tstep=5e-15;gpusteps=21;gpusave=20*tstep;runtime=gpusave;gpurun_number=1;savetstep=1;
  pulse_start=0;pulse_end=Inf;eval(src);
  mask=repmat(atomtype_==1,1,1,21);empty=~mask;norms=sqrt(mmx.^2+mmy.^2+mmz.^2);
  assert(max(abs(norms(mask)-1))<1e-12 && all(isfinite(norms),'all'));
  assert(all(mmx(empty)==0 & mmy(empty)==0 & mmz(empty)==0));
  assert(actual_updates==20 && numel(t)==21 && max(abs(t-(0:20)*tstep))<1e-25);
  trajectory=zeros(21,3*N);
  for j=1:21
   x=mmx(:,:,j);y=mmy(:,:,j);zspin=mmz(:,:,j);trajectory(j,:)=reshape([x(occ) y(occ) zspin(occ)]',1,[]);
  end
  f=fopen(fullfile(out,sprintf('%s-gpu%d.bin',tag,trial)),'w','ieee-le');assert(f>=0);
  fwrite(f,trajectory','double');fclose(f);
 end
 writejson(fullfile(out,[tag '-gpu-checks.json']),struct('field_error_T',ferr,'N',N,'updates',20,'vacancies_exact_zero',true));
 fprintf('Candidate source fields and GPU fixtures passed: %s, field error %.3g T\n',tag,ferr);
end
end
function v=extract(s,a,b)
i=strfind(s,a);j=strfind(s,b);assert(numel(i)==1);j=j(find(j>i,1));v=s(i:j-1);
end
function writejson(path,data)
f=fopen(path,'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(data,PrettyPrint=true));fclose(f);
end
