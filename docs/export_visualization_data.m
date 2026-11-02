function export_visualization_data(outputFile)
% CPU-only source snapshot and independent impulse checks. No integration.
% From the repository root:
% addpath('docs'); export_visualization_data('docs/mn3sn-data.js')
assert(nargin == 1, 'Specify the output JS path explicitly.');
root = fileparts(fileparts(mfilename('fullpath')));
original = pwd;
cleanup = onCleanup(@() cd(original));
cd(root);
names = {'main.m','constantfile.m','mn3sn_protocol.m','distrib.m','systemgeneration.m','integrate_llg.m'};
sources = struct(); hashes = struct();
for f = 1:numel(names)
    key = strrep(names{f},'.m','');
    sources.(key) = fileread(names{f});
    fid = fopen(names{f},'rb'); bytes = fread(fid,Inf,'*uint8'); fclose(fid);
    digest = java.security.MessageDigest.getInstance('SHA-256');
    digest.update(bytes);
    hashes.(key) = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
constantfile;
protocol=mn3sn_protocol('relaxation');
initial_theta_deg=protocol.seed_deg;
assert(strcmp(assignment(sources.main,'DMIenable'),'1'), ...
    'Exporter currently supports DMIenable=1 only; review the snapshot scope before refreshing.');
Jintra = assignment(sources.main,'Jintra');
Ksim = assignment(sources.main,'Ksim');
Dsim = assignment(sources.main,'Dsim');
% These expressions use ele from constantfile, never execute main.m.
Jintra = eval(Jintra); Ksim = eval(Ksim); Dsim = eval(Dsim);
assert(Jintra>0 && Ksim>0 && Dsim>0, 'Positive nonzero coefficients required by normalized impulse checks.');
musRE = 3*mub;
assert(strcmp(assignment(sources.main,'musRE'),'3*mub'), 'Moment mapping changed.');
text = sources.integrate_llg;
setup = between(text,'AsimnextL(2:2:natomW,:)=Jintra;','for ctL=1:natomL');
shifts = between(text,'        mmxnextL=circshift','        hex_x=');
fields = between(text,'        hex_x=','        if thermalenable');
qnames = {'nextL','previousL','nextW','previousW','previousWleft', ...
    'previousWright','nextWleft','nextWright'};
offsets = zeros(8,2);
for q = 1:8
    token = regexp(shifts,['mmx' qnames{q} '=circshift\(mmxtmp,\[(-?\d+),(-?\d+)\]\)'],'tokens','once');
    assert(~isempty(token), 'Shift expression changed.');
    offsets(q,:) = -[str2double(token{1}),str2double(token{2})];
end
% Independently transcribed coefficients of sz in By. Bz has opposite signs.
signs = [0 0 1 -1 0 1 0 -1; 1 1 0 -1 0 0 -1 0; -1 -1 1 0 1 0 0 0];
sizes = [5 5];
for W = [5 9 13]
    for L = [6 10 14], sizes(end+1,:) = [W L]; end %#ok<AGROW>
end
snapshots = cell(1,size(sizes,1));
for k = 1:size(sizes,1)
    natomW = sizes(k,1); natomL = sizes(k,2);
    load_fixed_atom_distrib = 0; save_fixed_atom_distrib = 0; dwcalc = 0;
    systemgeneration;
    for q = 1:8, eval(['Asim' qnames{q} '=zeros(natomW,natomL);']); end
    eval(setup);
    occupied = find(atomtype_==1);
    n = natomW*natomL;
    ex = zeros(n,8);
    for q = 1:8
        weights = eval(['Asim' qnames{q}]);
        ex(:,q) = weights(:)/Jintra;
    end
    dm = [hdmi_constant1(:) hdmi_constant2(:) hdmi_constant3(:)]*signs;
    axes = [coeffs_1(:) coeffs_2(:) coeffs_3(:)];
    spins = [mx_init(:) my_init(:) mz_init(:)];
    assert(all(all(spins(atomtype_(:)==3,:)==0)));
    assert(max(abs(sum(spins(occupied,:).^2,2)-1)) < 1e-12);
    assert(max(abs(sum(axes(occupied,:).^2,2)-1)) < 1e-12);
    if natomW==5 && natomL==5, assert(numel(occupied)==16); end
    if natomW==5 && natomL==6
        expected = [1 3 1 3 1 3;1 1 1 1 1 1;3 1 3 1 3 3;1 1 1 1 1 1;1 3 1 3 1 3];
        assert(isequal(atomtype_,expected));
    end
    s = struct('W',natomW,'L',natomL,'occupancy',atomtype_(:).', ...
        'axes',axes,'spins',spins,'scale',hex_constant(:).', ...
        'dmiMasks',[hdmi_constant1(:) hdmi_constant2(:) hdmi_constant3(:)], ...
        'exchange',ex,'dmi',dm);
    probes = cell(1,2);
    for bc = [1 0]
        % Build predicted operators from individual displayed reads.
        E = zeros(n); D = zeros(n);
        for i = occupied.'
            [r,c] = ind2sub([natomW natomL],i);
            for q = 1:8
                rr = r+offsets(q,1); cc = c+offsets(q,2);
                outside = rr<1 || rr>natomW || cc<1 || cc>natomL;
                if bc && outside, continue; end
                j = sub2ind([natomW natomL],mod(rr-1,natomW)+1,mod(cc-1,natomL)+1);
                if atomtype_(j)~=1, continue; end
                E(i,j) = E(i,j)+hex_constant(i)*ex(i,q);
                D(i,j) = D(i,j)+hdmi_constant(i)*dm(i,q);
            end
        end
        % Evaluate unmodified source expressions on each Cartesian impulse.
        measuredE = zeros(n); measuredD = zeros(n);
        muigpu = ones(natomW,natomL); muigpu(occupied)=musRE;
        for j = occupied.'
            for component = 1:3
                mmxtmp=zeros(natomW,natomL); mmytmp=mmxtmp; mmztmp=mmxtmp;
                if component==1, mmxtmp(j)=1;
                elseif component==2, mmytmp(j)=1;
                else, mmztmp(j)=1; end
                eval(shifts); eval(fields);
                efields = [hex_x(:) hex_y(:) hex_z(:)]*musRE/Jintra;
                target = zeros(n,3); target(:,component)=-E(:,j);
                closeEnough(efields(occupied,:),target(occupied,:));
                target=zeros(n,3);
                if component==2, target(:,3)=-D(:,j);
                elseif component==3, target(:,2)=D(:,j); end
                dfields=[zeros(n,1) hdmi_y(:) hdmi_z(:)]*musRE/Dsim;
                closeEnough(dfields(occupied,:),target(occupied,:));
                target=zeros(n,3); target(j,:)=axes(j,component)*axes(j,:);
                afields=[hani_x(:) hani_y(:) hani_z(:)]*musRE/(2*Ksim);
                closeEnough(afields(occupied,:),target(occupied,:));
                assert(hdmi_x==0);
                if component==1, measuredE(occupied,j)=-efields(occupied,1); end
                if component==3, measuredD(occupied,j)=dfields(occupied,2); end
            end
        end
        [ii,jj]=find(measuredE~=0 | measuredD~=0);
        idx=sub2ind([n n],ii,jj);
        mismatch=nnz(abs(E-E.')>1e-12)/2;
        edges=nnz(triu(E~=0 | E.'~=0,1));
        if natomW==5 && natomL==6 && bc==1
            assert(edges==30 && mismatch==12);
        end
        probes{2-bc}=struct('bc',bc,'operator',[ii-1 jj-1 measuredE(idx) measuredD(idx)], ...
            'exchangePairs',edges,'unequalPairs',mismatch,'maxExchangeAsymmetry',max(abs(E-E.'),[],'all'), ...
            'dmiAntisymmetryError',max(abs(D+D.'),[],'all'));
    end
    s.probes=probes;
    snapshots{k}=s;
    fprintf('%dx%d: %d occupied; both boundaries, exchange/DMI/anisotropy impulse checks passed\n',natomW,natomL,numel(occupied));
end
data=struct('schema',1,'description','Current-code snapshot; schematic coordinates only', ...
    'matlab',version,'hashes',hashes,'shiftNames',{qnames},'offsets',offsets, ...
    'parameters',struct('J',Jintra,'D',Dsim,'K',Ksim,'mu',musRE,'joulesPerMeV',ele*1e-3), ...
    'sources',sources,'snapshots',{snapshots},'tolerance',1e-12);
fid=fopen(outputFile,'w','n','UTF-8'); assert(fid~=-1,'Cannot open output.');
closer=onCleanup(@() fclose(fid));
fprintf(fid,'/* Generated by docs/export_visualization_data.m; do not hand edit. */\nwindow.MN3SN_DATA = %s;\n',jsonencode(data));
fprintf('Wrote verified static snapshot: %s\n',outputFile);
end

function value=assignment(text,name)
token=regexp(text,['(?m)^\s*' name '=([^;]+);'],'tokens','once');
assert(~isempty(token),'Parameter expression changed.'); value=token{1};
end

function part=between(text,first,last)
a=strfind(text,first); assert(numel(a)==1,'Ambiguous start marker.');
b=strfind(text(a:end),last); assert(~isempty(b),'Missing end marker.');
part=text(a:a+b(1)-2);
end

function closeEnough(a,b)
assert(all(isfinite(a),'all'),'Nonfinite static response.');
assert(max(abs(a-b),[],'all')<1e-12,'Viewer coefficient differs from source impulse response.');
assert(all(a(b==0)==0),'Unexpected nonzero response.');
end
