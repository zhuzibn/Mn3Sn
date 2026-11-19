function metadata=run_baseline_benchmark(mode,output_root,input_root)
%RUN_BASELINE_BENCHMARK Run a fresh deterministic suite; never overwrite output.
% mode: quick (4x4) or current (8x8). output_root must not exist.
% input_root optionally selects a prior suite supplying frozen inputs.
% Without it, require the accepted baseline_manifest.json; never bootstrap implicitly.
mode=validatestring(mode,{'quick','current'});
root=fileparts(fileparts(mfilename('fullpath')));
assert(nargin>=2 && ~isempty(output_root),'A fresh external output directory is required.');
bootstrap=nargin>=3 && isempty(input_root);
if nargin<3
    manifest_file=fullfile(root,'benchmark','baseline_manifest.json');
    assert(isfile(manifest_file),'Mn3Sn:MissingBaseline','Accepted baseline manifest is missing.');
    manifest=jsondecode(fileread(manifest_file));
    input_root=fullfile(root,manifest.(mode).directory);
end
if ~isempty(input_root),input_root=char(java.io.File(input_root).getCanonicalPath());end
frozen=cell(1,2);cases={'relaxation','driven_nonuniform'};
for k=1:2
    if ~isempty(input_root)
        data=load(fullfile(input_root,[cases{k} '.mat']),'result');
        benchmark_validate(data.result);frozen{k}=data.result;
    end
end
output_root=mn3sn_prepare_output(output_root,input_root,bootstrap);
original_path=path;original_dir=pwd;
cleanup=onCleanup(@() restore_context(original_dir,original_path));
diary(fullfile(output_root,'matlab.log'));diary_cleanup=onCleanup(@() diary('off'));
metadata=struct('schema_version',1,'mode',mode,'status','running', ...
    'timestamp',char(datetime('now','TimeZone','local','Format','yyyy-MM-dd HH:mm:ss Z')), ...
    'input_root',input_root,'matlab',version,'toolboxes',ver,'os',system_dependent('getos'));
metadata.command=sprintf('addpath(''benchmark''); run_baseline_benchmark(''%s'',''%s'',''%s'');', ...
    mode,strrep(output_root,'''',''''''),strrep(input_root,'''',''''''));
metadata.git_commit=git_output(root,'rev-parse HEAD');
metadata.git_status=git_output(root,'status --short --branch');
metadata.source=struct;source_root=fullfile(output_root,'source');mkdir(source_root);
files=[dir(fullfile(root,'*.m'));dir(fullfile(root,'benchmark','*.m'))];
for k=1:numel(files)
    file=fullfile(files(k).folder,files(k).name);
    if strcmp(files(k).folder,fullfile(root,'benchmark'))
        destination=fullfile(source_root,'benchmark');if ~isfolder(destination),mkdir(destination);end
        label=['benchmark_' files(k).name];
    else,destination=source_root;label=files(k).name;end
    copyfile(file,fullfile(destination,files(k).name));
    metadata.source.(matlab.lang.makeValidName(label))=file_hash(file);
end
% Production scripts execute exclusively in the isolated snapshot directory.
cd(source_root);addpath(source_root,'-begin');
try
    dev=gpuDevice(1);
    metadata.gpu=struct('name',dev.Name,'compute_capability',dev.ComputeCapability, ...
        'driver',dev.GraphicsDriverVersion,'toolkit',dev.ToolkitVersion, ...
        'total_memory',dev.TotalMemory);
    write_metadata(output_root,metadata);
    for k=1:2
        name=cases{k};fprintf('CASE %s / %s\n',mode,name);
        result=benchmark_case(mode,name,source_root,true,frozen{k});
        first_seconds=result.integration_seconds;warm=zeros(1,3);setup=zeros(1,4);
        fprintf('  first: %.6f s, norm error %.3g\n',first_seconds,result.observables.norm_error);
        setup(1)=result.setup_seconds;
        for repeat=1:3
            fresh=benchmark_case(mode,name,source_root,true,result);
            warm(repeat)=fresh.integration_seconds;setup(repeat+1)=fresh.setup_seconds;
            delta=max([max(abs(fresh.mmx-result.mmx),[],'all'), ...
                max(abs(fresh.mmy-result.mmy),[],'all'),max(abs(fresh.mmz-result.mmz),[],'all')]);
            assert(delta<=1e-10,'Timed repeat changed trajectory.');
            fprintf('  warm %d: %.6f s, max trajectory delta %.3g\n',repeat,warm(repeat),delta);
        end
        result.timing=struct('first_seconds',first_seconds,'warm_seconds',warm, ...
            'median_seconds',median(warm),'range_seconds',[min(warm) max(warm)], ...
            'setup_seconds',setup);
        result.provenance=metadata;
        io_timer=tic;save(fullfile(output_root,[name '.mat']),'result','-v7');
        plot_case(result,output_root);
        metadata.cases.(name)=result.timing;
        metadata.cases.(name).save_and_plot_seconds=toc(io_timer);
        metadata.cases.(name).result_sha256=file_hash(fullfile(output_root,[name '.mat']));
        metadata.cases.(name).norm_error=result.observables.norm_error;
        write_metadata(output_root,metadata);
    end
    metadata.status='complete';write_metadata(output_root,metadata);
    fprintf('PASS suite %s: %s\n',mode,output_root);
catch err
    metadata.status='failed';metadata.error=err.message;
    write_metadata(output_root,metadata);rethrow(err);
end
end

function plot_case(r,out)
f=figure('Visible','off','Color','w','Position',[100 100 1100 650]);
c=onCleanup(@() close(f));tiledlayout(2,3);
nexttile;imagesc(r.atomtype);axis image;colorbar;title('Occupation mask (array indices)');
nexttile;imagesc(r.mmz(:,:,1));axis image;colorbar;title('Initial m_z');
nexttile;imagesc(r.mmz(:,:,end));axis image;colorbar;title('Final m_z');
nexttile;plot(r.t/1e-12,r.observables.mean);xlabel('Time (ps)');ylabel('Occupied mean');legend('x','y','z');
nexttile;plot(r.t/1e-12,squeeze(r.observables.class_mean(:,3,:)));xlabel('Time (ps)');ylabel('Class mean z');legend('a','b','c');
nexttile;n=sqrt(r.mmx.^2+r.mmy.^2+r.mmz.^2);mask=repmat(r.atomtype==1,1,1,size(n,3));
histogram(n(mask)-1);xlabel('Occupied norm error (all saved states)');
sgtitle([r.mode ' / ' strrep(r.case_name,'_',' ')]);
exportgraphics(f,fullfile(out,[r.case_name '.png']),'Resolution',140);
end

function value=git_output(root,args)
[status,value]=system(sprintf('git -C "%s" %s',root,args));
assert(status==0,'Git provenance failed: %s',value);value=strtrim(value);
end

function value=file_hash(file)
f=fopen(file,'rb');assert(f>=0);c=onCleanup(@() fclose(f));bytes=fread(f,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function write_metadata(out,data)
f=fopen(fullfile(out,'metadata.json'),'w');assert(f>=0);c=onCleanup(@() fclose(f));
fprintf(f,'%s\n',jsonencode(data,PrettyPrint=true));
end

function restore_context(original_dir,original_path)
cd(original_dir);path(original_path);
end
