function check_three_spin(out)
% Source-based static gates; no production main, output, or GPU integration.
original=pwd; cleanup=onCleanup(@() cd(original));
addpath(fullfile(out,'native')); addpath(fullfile(out,'reference'));
cd(out);
native=parameters(fullfile(out,'native'));
cd(fullfile(out,'reference')); initialize;
ref=struct('mu',3*mub,'J',Hex_2*3*mub,'D',DMI*3*mub, ...
    'Kfield',hk6*3*mub,'gamma',gam,'alpha',alp,'Ms',Ms_Fe*1e3, ...
    'thickness',tFL,'b',thetaSH*5e10/(2*tFL*Ms_Fe*1e3/hbar));
cd(out);
axes=[zeros(3,1),sin([1;9;5]*pi/6),cos([1;9;5]*pi/6)];
% Stage 1: original occupation, initializer, and all 9 Cartesian impulses.
[~,~,~,geometry]=spatial(zeros(3),native,fullfile(out,'native'));
assert(isequal(geometry.occupation,[1 3;1 1]));
assert(max(abs(geometry.axes-axes),[],'all')<1e-14);
assert(all(geometry.initial(3,:)==0)); % column-major empty slot (1,2)
assert(max(abs(vecnorm(geometry.initial([4 1 2],:),2,2)-1))<1e-12);
operators=zeros(9,9,3,2); fieldError=zeros(1,3);
for parameter=1:2
    if parameter==1, p=native; else, p=ref; end
    for q=1:9
        m=zeros(3); m(q)=1;
        [ex,dm,an]=spatial(m,p,fullfile(out,'native'));
        operators(:,q,1,parameter)=ex(:);
        operators(:,q,2,parameter)=dm(:);
        operators(:,q,3,parameter)=an(:);
        expected=triangle(m,p,axes);
        measured=cat(3,ex,dm,an);
        for term=1:3
            fieldError(term)=max(fieldError(term),max(abs(measured(:,:,term)-expected(:,:,term)),[],'all'));
        end
    end
end
E=-operators(1:3,1:3,1,1)*native.mu/native.J;
Y=operators(4:6,7:9,2,1)*native.mu/native.D;
assert(max(abs(E-(ones(3)-eye(3))),[],'all')<1e-14);
assert(max(abs(Y+Y'),[],'all')<1e-14);
assert(isequal(round(Y),[0 -1 1;1 0 -1;-1 1 0]));
assert(all(fieldError<1e-11));
fprintf('PASS mapping/topology: 3 occupied, reciprocal exchange, cyclic antisymmetric DMI\n');
% Stage 2: reference field_eta term isolation and independent energy gradient.
rng(20261102); fdError=zeros(2,4); referenceError=zeros(1,4);
for trial=1:8
    m=randn(3); m=m./vecnorm(m,2,2);
    B=[.013 .025 -.019];
    expected=triangle(m,ref,axes);
    for term=1:4
        termp=ref; termp.J=ref.J*(term==1); termp.D=ref.D*(term==2); termp.Kfield=ref.Kfield*(term==3);
        external=B*(term==4);
        h=reference_fields(m,termp,external,0);
        if term<=3, target=expected(:,:,term); else, target=repmat(B,3,1); end
        referenceError(term)=max(referenceError(term),max(abs(h-target),[],'all'));
    end
    for parameter=1:2
        if parameter==1, p=native; else, p=ref; end
        target=cat(3,triangle(m,p,axes),repmat(B,3,1));
        for term=1:4
            fd=zeros(3); delta=1e-5;
            for q=1:9
                a=m; b=m; a(q)=a(q)+delta; b(q)=b(q)-delta;
                fd(q)=-(energy(a,p,axes,B,term)-energy(b,p,axes,B,term))/(2*delta);
            end
            % Energy divided by mu gives a well-scaled tesla-valued gradient.
            residual=max(abs(fd-target(:,:,term)),[],'all')/max(1,max(abs(target(:,:,term)),[],'all'));
            fdError(parameter,term)=max(fdError(parameter,term),residual);
        end
    end
end
assert(all(referenceError<1e-11)); assert(all(fdError<1e-7,'all'));
fprintf('PASS fields/energy: source discrepancy %.3g T; scaled gradient residual %.3g\n',max(referenceError),max(fdError,[],'all'));
save(fullfile(out,'static.mat'),'native','ref','axes','operators','geometry');
result=struct('pass',true,'matlab',version,'toolboxes',{ver},'native',native,'reference',ref, ...
    'geometry',geometry,'exchange',E,'dmi_By_from_z',Y,'field_error_T',fieldError, ...
    'reference_field_error_T',referenceError,'gradient_relative_error',fdError);
writejson(fullfile(out,'static.json'),result);
end

function p=parameters(folder)
original=pwd; cleanup=onCleanup(@() cd(original)); cd(folder);
constantfile;
src=fileread('main.m');
names={'Jintra','Ksim','Dsim','gRE','volume','d','alp','thetaSH'};
for i=1:numel(names)
    token=regexp(src,['(?m)^\s*' names{i} '=([^;]+);'],'tokens','once');
    % gRE follows gTM on the same source line.
    if isempty(token), token=regexp(src,[';' names{i} '=([^;]+);'],'tokens','once'); end
    assert(~isempty(token),'Missing source assignment: %s',names{i}); eval([names{i} '=' token{1} ';']);
end
p=struct('mu',3*mub,'J',Jintra,'D',Dsim,'Kfield',Ksim, ...
    'gamma',gRE*mub/(hbar*ele),'alpha',alp,'Ms',3*mub/volume, ...
    'thickness',d,'b',hbar/2*thetaSH*5e10/(3*mub/volume*d));
end

function [ex,dm,an,geometry]=spatial(m,p,folder)
% Execute original setup/stencils, not a transcribed triangle field.
natomW=2; natomL=2; bc=1; dwcalc=0;
load_fixed_atom_distrib=0; save_fixed_atom_distrib=0;
eval(fileread(fullfile(folder,'distrib.m')));
% systemgeneration calls distrib from the isolated native source path.
eval(fileread(fullfile(folder,'systemgeneration.m')));
map=[4 1 2];
geometry=struct('occupation',atomtype_,'axes',[coeffs_1(map)' coeffs_2(map)' coeffs_3(map)'], ...
    'initial',[mx_init(:) my_init(:) mz_init(:)],'map',map);
Jintra=p.J; Ksim=p.Kfield; Dsim=p.D;
muigpu=ones(2); muigpu(atomtype_==1)=p.mu;
src=fileread(fullfile(folder,'integrate_llg.m'));
setup=between(src,'AsimnextL=zeros','for ctL=1:natomL');
eval(setup);
mmxtmp=zeros(2); mmytmp=zeros(2); mmztmp=zeros(2);
mmxtmp(map)=m(:,1); mmytmp(map)=m(:,2); mmztmp(map)=m(:,3);
eval(between(src,'        mmxnextL=circshift','        if thermalenable'));
ex=[hex_x(map)' hex_y(map)' hex_z(map)'];
dm=[zeros(3,1) hdmi_y(map)' hdmi_z(map)'];
an=[hani_x(map)' hani_y(map)' hani_z(map)'];
end

function terms=triangle(m,p,e)
ex=-p.J/p.mu*(m([2 1 1],:)+m([3 3 2],:));
dm=p.D/p.mu*cross(repmat([1 0 0],3,1),m([2 3 1],:)-m([3 1 2],:),2);
an=p.Kfield/p.mu*sum(m.*e,2).*e;
terms=cat(3,ex,dm,an);
end

function u=energy(m,p,e,B,term)
if term==1, u=p.J/p.mu*sum(sum(m.*m([2 3 1],:)));
elseif term==2
    c=cross(m,m([2 3 1],:),2); u=p.D/p.mu*sum(c(:,1));
elseif term==3, u=-p.Kfield/(2*p.mu)*sum(sum(m.*e,2).^2);
else, u=-sum(m*B'); end
end

function [h,b]=reference_fields(m,p,B,current)
[a,bv,c,~,~,b,~,~,~,~,~,~,~,~]=field_eta(p.Kfield/p.mu,m(1,:),m(2,:),m(3,:), ...
    [0 0 1],[0 0 1],[0 0 1],0,0,0,0,zeros(3),B,0, ...
    p.thickness,p.Ms/1e3,p.Ms/1e3,0,zeros(3),[0 0 1],0,50e-9,50e-9,0, ...
    .06,7e-9,5e-9,current,0,p.alpha,5e-15,0,0,p.J/p.mu,0,0,0,p.D/p.mu);
h=[a;bv;c];
end

function part=between(src,first,last)
a=strfind(src,first); assert(numel(a)==1);
b=strfind(src(a:end),last); assert(~isempty(b)); part=src(a:a+b(1)-2);
end

function writejson(path,data)
f=fopen(path,'w'); assert(f>=0); cleanup=onCleanup(@() fclose(f));
fprintf(f,'%s\n',jsonencode(data,PrettyPrint=true));
end
