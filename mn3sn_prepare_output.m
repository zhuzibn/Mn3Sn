function output_root=mn3sn_prepare_output(output_root,input_root,bootstrap)
%MN3SN_PREPARE_OUTPUT Claim a fresh external directory; never replace results.
% input_root additionally protects an explicit benchmark input tree. bootstrap
% is reserved for separately authorized reference generation, never ordinary runs.
if nargin<2,input_root='';end
if nargin<3,bootstrap=false;end
root=fileparts(mfilename('fullpath'));
assert((ischar(output_root) && isrow(output_root)) || ...
    (isstring(output_root) && isscalar(output_root) && ~ismissing(output_root)), ...
    'Mn3Sn:OutputPath','Output must be a nonempty path to a fresh directory.');
output_root=char(output_root);
assert(~isempty(strtrim(output_root)),'Mn3Sn:OutputPath','Output path must not be empty.');
output_root=canonical(output_root);root=canonical(root);
assert(~within(output_root,root),'Mn3Sn:OutputPath','Outputs must be outside the source repository.');
assert(~contains(lower(output_root),'onedrive'),'Mn3Sn:OutputPath','Outputs must be outside OneDrive.');
manifest_file=fullfile(root,'benchmark','baseline_manifest.json');
assert(isfile(manifest_file) || bootstrap,'Mn3Sn:MissingBaseline', ...
    'Missing baseline_manifest.json; ordinary writes cannot establish baseline protection.');
protected={};
if isfile(manifest_file)
    manifest=jsondecode(fileread(manifest_file));
    assert(manifest.schema_version==1,'Mn3Sn:MissingBaseline','Unsupported baseline manifest.');
    for mode={'quick','current'}
        folder=manifest.(mode{1}).directory;
        assert(ischar(folder) && isrow(folder) && ~isempty(folder), ...
            'Mn3Sn:MissingBaseline','Invalid registered baseline directory.');
        protected{end+1}=canonical(fullfile(root,folder)); %#ok<AGROW>
    end
end
if ~isempty(input_root),protected{end+1}=canonical(char(input_root));end
for k=1:numel(protected)
    assert(~within(output_root,protected{k}),'Mn3Sn:ProtectedBaseline', ...
        'Output is inside a protected baseline: %s.',protected{k});
end
assert(~isfolder(output_root) && ~isfile(output_root),'Mn3Sn:OutputExists', ...
    'Output already exists: %s.',output_root);
parent=fileparts(output_root);if ~isfolder(parent),mkdir(parent);end
% mkdir() itself is the exclusive claim, including competing directory creation.
assert(java.io.File(output_root).mkdir(),'Mn3Sn:OutputExists', ...
    'Could not claim a fresh output directory: %s.',output_root);
end

function value=canonical(value)
value=char(java.io.File(value).getCanonicalPath());
end

function yes=within(child,parent)
yes=strcmpi(child,parent) || startsWith(lower(child),[lower(parent) filesep]);
end
