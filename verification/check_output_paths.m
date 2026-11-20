function check_output_paths(out)
% Real filesystem aliases in disposable fixtures; no accepted tree is written.
root=fileparts(fileparts(mfilename('fullpath')));
assert(~isfolder(out));mkdir(out);
oldcwd=pwd;oldpath=path;cleanup=onCleanup(@() restore(oldcwd,oldpath));
source=fullfile(out,'source');mkdir(fullfile(source,'benchmark'));
copyfile(fullfile(root,'mn3sn_prepare_output.m'),source);
registered=fullfile(out,'registered');selected=fullfile(out,'selected');
cloud=fullfile(out,'OneDrive-fixture');mkdir(registered);mkdir(selected);mkdir(cloud);
manifest=struct('schema_version',1,'quick',struct('directory',fullfile('..','registered')),...
    'current',struct('directory',fullfile('..','registered')));
f=fopen(fullfile(source,'benchmark','baseline_manifest.json'),'w');
assert(f>=0);fprintf(f,'%s',jsonencode(manifest));fclose(f);
cd(source);addpath(source,'-begin');clear mn3sn_prepare_output;rehash;
paths={fullfile(source,'fresh'),fullfile('benchmark','fresh'),...
    fullfile(registered,'nested','fresh'),fullfile(selected,'fresh'),fullfile(cloud,'fresh')};
ids={'Mn3Sn:OutputPath','Mn3Sn:OutputPath','Mn3Sn:ProtectedBaseline',...
    'Mn3Sn:ProtectedBaseline','Mn3Sn:OutputPath'};
if ispc
    targets={source,registered,selected,cloud};
    for k=1:numel(targets)
        alias=fullfile(out,sprintf('junction-%d',k));
        [status,message]=system(sprintf('cmd /c mklink /J "%s" "%s"',alias,targets{k}));
        assert(status==0,'Mn3Sn:TestFailure','Cannot create junction: %s',message);
        paths{end+1}=fullfile(alias,'nested','fresh'); %#ok<AGROW>
    end
    ids=[ids,{'Mn3Sn:OutputPath','Mn3Sn:ProtectedBaseline',...
        'Mn3Sn:ProtectedBaseline','Mn3Sn:OutputPath'}];
    paths{end+1}=upper(fullfile(registered,'fresh'));ids{end+1}='Mn3Sn:ProtectedBaseline';
end
for k=1:numel(paths)
    try
        mn3sn_prepare_output(paths{k},selected);
        error('Mn3Sn:TestFailure','Protected destination accepted.');
    catch e
        assert(strcmp(e.identifier,ids{k}),'Mn3Sn:TestFailure','Wrong rejection: %s',e.message);
    end
    assert(~isfolder(paths{k}));
end
fresh=fullfile(out,'ordinary','nested','fresh');
assert(isfolder(mn3sn_prepare_output(fresh,selected)));
try
    mn3sn_prepare_output(fresh,selected);error('Mn3Sn:TestFailure','Existing output accepted.');
catch e,assert(strcmp(e.identifier,'Mn3Sn:OutputExists'));end
fprintf('PASS output aliases: %d protected paths, fresh nested output and existing rejection\n',numel(paths));
end

function restore(folder,oldpath)
cd(folder);path(oldpath);clear mn3sn_prepare_output;rehash;
end
