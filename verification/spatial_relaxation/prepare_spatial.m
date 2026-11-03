function prepare_spatial(out)
root=pwd; addpath(root); dev=gpuDevice(1); constantfile;
src=fileread('integrate_llg.m');
setup=extract(src,'AsimnextL(2:2:natomW,:)=Jintra;','for ctL=1:natomL');
shifts=extract(src,'        mmxnextL=circshift','        hex_x=');
fields=extract(src,'        hex_x=','        if thermalenable');
sizes=[2 2;5 5;5 6;9 10;13 14;21 22];
meta=struct('matlab',version,'gpu',dev.Name,'sizes',sizes); writejson(fullfile(out,'environment.json'),meta);
for q=1:size(sizes,1)
 natomW=sizes(q,1);natomL=sizes(q,2);natom=natomW*natomL;
 dwcalc=0;load_fixed_atom_distrib=0;save_fixed_atom_distrib=0;initial_theta_deg=121;
 systemgeneration; occ=find(atomtype_==1);N=numel(occ);bc=1;
 Jintra=100.812*(3*mub);Ksim=.196e-3*ele;Dsim=.833e-3*ele;musRE=3*mub;gamRE=gam;alp=.003;
 names={'nextL','previousL','nextW','previousW','previousWleft','previousWright','nextWleft','nextWright'};
 for k=1:8,eval(['Asim' names{k} '=zeros(natomW,natomL);']);end
 eval(setup);muigpu=ones(natomW,natomL);muigpu(occ)=musRE;
 A=zeros(N);C=zeros(N);axes=[coeffs_1(occ) coeffs_2(occ) coeffs_3(occ)];
 seed=[mx_init(occ) my_init(occ) mz_init(occ)];
 for j=1:N
  mmxtmp=zeros(natomW,natomL);mmytmp=mmxtmp;mmztmp=mmxtmp;
  mmxtmp(occ(j))=1;mmztmp(occ(j))=1;eval(shifts);eval(fields);
  A(:,j)=-hex_x(occ)*musRE/Jintra;C(:,j)=hdmi_y(occ)*musRE/Dsim;
 end
 % Independent random-state field test against the linear operator.
 rng(20261102);r=randn(N,3);r=r./vecnorm(r,2,2);
 mmxtmp=zeros(natomW,natomL);mmytmp=mmxtmp;mmztmp=mmxtmp;
 mmxtmp(occ)=r(:,1);mmytmp(occ)=r(:,2);mmztmp(occ)=r(:,3);eval(shifts);eval(fields);
 predicted=-Jintra/musRE*A*r+Dsim/musRE*[zeros(N,1) C*r(:,3) -C*r(:,2)]+2*Ksim/musRE*sum(r.*axes,2).*axes;
 actual=[hex_x(occ)+hani_x(occ),hex_y(occ)+hani_y(occ)+hdmi_y(occ),hex_z(occ)+hani_z(occ)+hdmi_z(occ)];
 ferr=max(abs(predicted-actual),[],'all');assert(ferr<1e-11);
 tag=sprintf('%dx%d',natomW,natomL);
 [ii,jj]=find(A~=0|C~=0);idx=sub2ind([N N],ii,jj);
 f=fopen(fullfile(out,[tag '.model']),'w');assert(f>=0);
 fprintf(f,'%d %d %d %.17g %.17g %.17g %.17g %.17g\n',N,natomW,natomL,Jintra/musRE,Dsim/musRE,2*Ksim/musRE,alp,gamRE);
 fprintf(f,'%.17g %.17g %.17g %.17g %.17g %.17g\n',[axes seed]');
 fprintf(f,'%d\n',numel(ii));fprintf(f,'%d %d %.17g %.17g\n',[ii-1 jj-1 A(idx) C(idx)]');fclose(f);
 stats=struct('occupied',N,'pairs',nnz(triu(A~=0|A'~=0)), ...
  'exchange_asymmetry',max(abs(A-A'),[],'all'),'dmi_antisymmetry',max(abs(C+C'),[],'all'), ...
  'unequal_exchange_pairs',nnz(abs(A-A')>1e-12)/2,'field_error_T',ferr, ...
  'unweighted_exchange_reciprocal',isequal(A~=0,A'~=0), ...
  'unweighted_dmi_antisymmetric',isequal(sign(C),-sign(C')));
 writejson(fullfile(out,[tag '-static.json']),stats);
 save(fullfile(out,[tag '-static.mat']),'A','C','axes','seed','occ','atomtype_','hex_constant');
 % Actual production GPU: 20 updates from seed and a noncoplanar perturbation.
 original_scale=hex_constant;
 for trial=0:3
  hex_constant=original_scale;hdmi_constant=original_scale;
  if trial>=2,hex_constant(:)=.5;hdmi_constant(:)=.5;end
  m=seed;if mod(trial,2),m=m+.01*r;m=m./vecnorm(m,2,2);end
  mx_init(occ)=m(:,1);my_init(occ)=m(:,2);mz_init(occ)=m(:,3);
  BDSOTRE=0;BDSTTRE=0;chi=0;dipolemode=0;thermalenable=0;enablefixedge=0;rk4=1;
  psjSHEx=1;psjSHEy=0;psjSHEz=0;psjSTTx=0;psjSTTy=0;psjSTTz=1;
  Hext=[0 0 0];tstep=5e-15;gpusteps=21;gpusave=20*tstep;runtime=gpusave;gpurun_number=1;savetstep=1;
  pulse_start=0;pulse_end=Inf;eval(src);
  mask=repmat(atomtype_==1,1,1,21);empty=~mask;norms=sqrt(mmx.^2+mmy.^2+mmz.^2);
  assert(max(abs(norms(mask)-1))<1e-12 && all(isfinite(norms),'all'));
  assert(all(mmx(empty)==0 & mmy(empty)==0 & mmz(empty)==0));assert(actual_updates==20);
  trajectory=zeros(21,3*N);
  for j=1:21
   x=mmx(:,:,j);y=mmy(:,:,j);z=mmz(:,:,j);trajectory(j,:)=reshape([x(occ) y(occ) z(occ)]',1,[]);
  end
  f=fopen(fullfile(out,sprintf('%s-gpu%d.bin',tag,trial)),'w','ieee-le');fwrite(f,trajectory','double');fclose(f);
 end
 fprintf('Prepared %s: %d spins, %d unequal pairs, field error %.3g T\n',tag,N,stats.unequal_exchange_pairs,ferr);
end
end
function v=extract(s,a,b)
i=strfind(s,a);j=strfind(s,b);assert(numel(i)==1);j=j(find(j>i,1));v=s(i:j-1);
end
function writejson(path,data)
f=fopen(path,'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(data,PrettyPrint=true));fclose(f);
end
