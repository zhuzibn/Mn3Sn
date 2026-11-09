function prepare_periodic(out)
% Export and verify the production periodic model, without source-text overrides.
root=pwd;addpath(root);constantfile;dev=gpuDevice(1);
assert(isfolder(out) && ~isfile(fullfile(out,'gpu-gates.json')),'Use a fresh output directory.');
src=strrep(fileread('integrate_llg.m'),sprintf('\r\n'),newline);
distribution=strrep(fileread('distrib.m'),sprintf('\r\n'),newline);
generation=strrep(fileread('systemgeneration.m'),sprintf('\r\n'),newline);
writefile(fullfile(out,'periodic-distrib.txt'),distribution);
writefile(fullfile(out,'periodic-integrator.txt'),src);
writejson(fullfile(out,'environment.json'),struct('matlab',version,'gpu',dev.Name, ...
 'compute_capability',dev.ComputeCapability,'toolboxes',{ver}));
sizes=[2 2;4 4;8 8;12 12;20 20];checks={};
for q=1:size(sizes,1)
 natomW=sizes(q,1);natomL=sizes(q,2);natom=natomW*natomL;
 dwcalc=0;load_fixed_atom_distrib=0;save_fixed_atom_distrib=0;initial_theta_deg=30;
 bc=double(q==1); % Open three-spin control; periodic cells use production directly.
 eval(generation);integration_source=src;
 occ=find(atomtype_==1);N=numel(occ);tag=sprintf('%dx%d',natomW,natomL);
 Jintra=100.812*(3*mub);Ksim=.196e-3*ele;Dsim=.833e-3*ele;musRE=3*mub;gamRE=gam;alp=.003;
 names={'nextL','previousL','nextW','previousW','previousWleft','previousWright','nextWleft','nextWright'};
 for k=1:8,eval(['Asim' names{k} '=zeros(natomW,natomL);']);end
 eval(extract(integration_source,'AsimnextL(2:2:natomW,:)=Jintra;','for ctL=1:natomL'));
 muigpu=ones(natomW,natomL);muigpu(occ)=musRE;
 field_source=extract(integration_source,'        mmxnextL=circshift','        if thermalenable');
 A=zeros(N);C=zeros(N);axes=[coeffs_1(occ) coeffs_2(occ) coeffs_3(occ)];
 seed=[mx_init(occ) my_init(occ) mz_init(occ)];
 for j=1:N
  mmxtmp=zeros(natomW,natomL);mmytmp=mmxtmp;mmztmp=mmxtmp;
  mmxtmp(occ(j))=1;mmztmp(occ(j))=1;eval(field_source);
  A(:,j)=-hex_x(occ)*musRE/Jintra;C(:,j)=hdmi_y(occ)*musRE/Dsim;
 end
 assert(max(abs(A-A'),[],'all')<1e-15 && max(abs(C+C'),[],'all')<1e-15);
 assert(isequal(A~=0,C~=0) && all(diag(A)==0));
 classes=ones(N,1);classes(axes(:,2)<-.5)=2;classes(axes(:,3)<0 & axes(:,2)>-.5)=3;
 counts=accumarray(classes,1,[3 1]);assert(all(counts==N/3));
 if q>1,assert(N==3*natom/4 && all(sum(A~=0,2)==4));end
 targetA=[0 1 1;1 0 1;1 1 0];targetC=[0 -1 1;1 0 -1;-1 1 0];
 for k=1:3
  assert(max(abs(sum(A(:,classes==k),2)-targetA(classes,k)))<1e-15);
  assert(max(abs(sum(C(:,classes==k),2)-targetC(classes,k)))<1e-15);
 end
 rng(20261109);r=randn(N,3);r=r./vecnorm(r,2,2);
 mmxtmp=zeros(natomW,natomL);mmytmp=mmxtmp;mmztmp=mmxtmp;
 mmxtmp(occ)=r(:,1);mmytmp(occ)=r(:,2);mmztmp(occ)=r(:,3);eval(field_source);
 predicted=-Jintra/musRE*A*r+Dsim/musRE*[zeros(N,1) C*r(:,3) -C*r(:,2)]+2*Ksim/musRE*sum(r.*axes,2).*axes;
 actual=[hex_x(occ)+hani_x(occ),hex_y(occ)+hani_y(occ)+hdmi_y(occ),hex_z(occ)+hani_z(occ)+hdmi_z(occ)];
 ferr=max(abs(predicted-actual),[],'all');assert(ferr<1e-11);
 derivative=zeros(N,3);
 for j=1:N
  [derivative(j,1),derivative(j,2),derivative(j,3)]=atom_llgs_rhs(r(j,1),r(j,2),r(j,3),1,0,0,0,0,1, ...
   gam/(1+alp^2),alp,actual(j,1),actual(j,2),actual(j,3),0,0,0,0);
 end
 residual=max(abs(derivative+gam*cross(r,actual,2)-alp*cross(r,derivative,2)),[],'all')/gam;
 assert(residual<1e-10);
 [ii,jj]=find(A~=0|C~=0);idx=sub2ind([N N],ii,jj);
 f=fopen(fullfile(out,[tag '.model']),'w');assert(f>=0);
 fprintf(f,'%d %d %d %.17g %.17g %.17g %.17g %.17g\n',N,natomW,natomL,Jintra/musRE,Dsim/musRE,2*Ksim/musRE,alp,gamRE);
 fprintf(f,'%.17g %.17g %.17g %.17g %.17g %.17g\n',[axes seed]');
 fprintf(f,'%d\n',numel(ii));fprintf(f,'%d %d %.17g %.17g\n',[ii-1 jj-1 A(idx) C(idx)]');fclose(f);
 [rows,cols]=ind2sub([natomW natomL],occ);
 save(fullfile(out,[tag '-static.mat']),'A','C','axes','seed','occ','atomtype_','classes','rows','cols');
 writejson(fullfile(out,[tag '-sites.json']),struct('W',natomW,'L',natomL,'rows',rows,'cols',cols,'classes',classes-1));
 for trial=0:1
  m=seed;if trial,m=m+.01*r;m=m./vecnorm(m,2,2);end
  mx_init(occ)=m(:,1);my_init(occ)=m(:,2);mz_init(occ)=m(:,3);
  BDSOTRE=0;BDSTTRE=0;chi=0;dipolemode=0;thermalenable=0;enablefixedge=0;rk4=1;
  psjSHEx=1;psjSHEy=0;psjSHEz=0;psjSTTx=0;psjSTTy=0;psjSTTz=1;
  Hext=[0 0 0];tstep=5e-15;gpusteps=11;gpusave=10*tstep;runtime=20*tstep;gpurun_number=2;savetstep=1;
  pulse_start=0;pulse_end=Inf;eval(integration_source);
  mask=repmat(atomtype_==1,1,1,21);empty=~mask;norms=sqrt(mmx.^2+mmy.^2+mmz.^2);
  assert(max(abs(norms(mask)-1))<1e-12 && all(isfinite(norms),'all'));
  assert(all(mmx(empty)==0 & mmy(empty)==0 & mmz(empty)==0));
  assert(actual_updates==20 && numel(t)==21 && max(abs(t-(0:20)*tstep))<1e-25);
  trajectory=zeros(21,3*N);
  for j=1:21
   x=mmx(:,:,j);y=mmy(:,:,j);z=mmz(:,:,j);trajectory(j,:)=reshape([x(occ) y(occ) z(occ)]',1,[]);
  end
  f=fopen(fullfile(out,sprintf('%s-gpu%d.bin',tag,trial)),'w','ieee-le');assert(f>=0);
  fwrite(f,trajectory','double');fclose(f);
 end
 checks{end+1}=struct('size',tag,'N',N,'field_error_T',ferr,'gilbert_residual_over_gamma_T',residual, ...
  'counts',counts,'updates',20,'vacancies_exact_zero',true); %#ok<AGROW>
 fprintf('Periodic preparation passed %s: %d spins, field error %.3g T\n',tag,N,ferr);
end
writejson(fullfile(out,'gpu-gates.json'),[checks{:}]);
end
function v=extract(s,a,b)
i=strfind(s,a);j=strfind(s,b);assert(numel(i)==1);j=j(find(j>i,1));assert(~isempty(j));v=s(i:j-1);
end
function writejson(path,data)
writefile(path,[jsonencode(data,PrettyPrint=true) newline]);
end
function writefile(path,data)
f=fopen(path,'w');assert(f>=0);fprintf(f,'%s',data);fclose(f);
end
