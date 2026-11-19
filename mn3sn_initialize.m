function state=mn3sn_initialize(cfg)
%MN3SN_INITIALIZE Pure CPU initialization and exact legacy coefficient mapping.
root=fileparts(mfilename('fullpath'));original=pwd;original_path=path;
cleanup=onCleanup(@() restore(original,original_path));
cd(root);addpath(root,'-begin');
cfg=mn3sn_validate_config(cfg);
constantfile;
% Only validated known inputs become legacy script aliases.
names=setdiff(fieldnames(cfg),{'protocol','initial_state','gpu_device','rng_seed'});
for k=1:numel(names),eval([names{k} '=cfg.(names{k});']);end
protocol=cfg.protocol;
natomW=protocol.natomW;natomL=protocol.natomL;initial_theta_deg=protocol.seed_deg;
natom=natomW*natomL;
systemgeneration;
Jfefe=-2.835e-21;Jfegd=1.09e-21;gTM=2.2;
gamTM=gTM*mub/(hbar*ele);musTM=2.217*mub;msTM=musTM/volume;
if ~DMIenable,Dsim=0;end
jcSOT=protocol.current;Hext=protocol.field;
pulse_start=protocol.pulse_start;pulse_end=protocol.pulse_end;
BDSOTRE=SOT_DLT*hbar/2*thetaSH*jcSOT/(msRE*tz);
BDSOTTM=SOT_DLT*hbar/2*thetaSH*jcSOT/(msTM*tz);
BDSTTRE=STT_DLT*hbar/2*etaSTT*jcSTT/(msRE*tz);
BDSTTTM=STT_DLT*hbar/2*etaSTT*jcSTT/(msTM*tz);
assert(isfinite(msRE*tz) && msRE*tz>0 && isfinite(msTM*tz) && msTM*tz>0, ...
    'Mn3Sn:Configuration','msRE*tz and derived msTM*tz must be finite and positive.');
assert(isfinite(1+alp^2) && isfinite(gamRE/(1+alp^2)) && gamRE/(1+alp^2)>0, ...
    'Mn3Sn:Configuration','alp and gamRE must give a finite positive Gilbert scaling.');
assert(all(isfinite([Jintra/musRE Dsim/musRE 2*Ksim/musRE BDSOTRE BDSOTTM BDSTTRE BDSTTTM])), ...
    'Mn3Sn:Configuration','Material, moment and drive inputs must give finite field/torque coefficients.');
gpusave=protocol.gpusave;runtime=protocol.runtime;tstep=protocol.tstep;savetstep=protocol.savetstep;
gpurun_number=round(runtime/gpusave);gpusteps=round(gpusave/tstep)+1;
totstep=gpurun_number*(gpusteps-1);
names=strsplit(['gam ele mu_0 mub hbar kb rk4 bc DMIenable dwcalc thermalenable ' ...
    'loadstartm load_fixed_atom_distrib save_fixed_atom_distrib dipolemode enablefixedge ' ...
    'natomW natomL initial_theta_deg d natom Ksim Jintra Jfefe Jfegd gTM gamTM gamRE ' ...
    'tz volume musRE musTM msRE msTM Dsim alp jcSOT jcSTT Hext pulse_start pulse_end ' ...
    'SOT_DLT SOT_FLT psjSHE thetaSH chi BDSOTRE BDSOTTM STT_DLT STT_FLT psjSTT ' ...
    'etaSTT chiSTT BDSTTRE BDSTTTM T gpusave runtime tstep savetstep gpurun_number gpusteps totstep']);
state.config=struct;
for k=1:numel(names),state.config.(names{k})=eval(names{k});end
state.config.protocol=protocol;
state.atomtype=atomtype_;state.axes=cat(3,coeffs_1,coeffs_2,coeffs_3);
state.generated_initial=cat(3,mx_init,my_init,mz_init);
state.initial=state.generated_initial;
if ~isempty(cfg.initial_state),state.initial=cfg.initial_state;end
occ=atomtype_==1;empty=~occ;x=state.initial;
assert(all(x(repmat(empty,1,1,3))==0),'Mn3Sn:InitialState','Initial vacant slots must be exactly zero.');
norms=sqrt(sum(x.^2,3));
assert(all(isfinite(norms(occ))) && max(abs(norms(occ)-1))<=1e-12, ...
    'Mn3Sn:InitialState','Initial occupied spins must have unit norm within 1e-12.');
[row,col]=ndgrid(1:natomW,1:natomL);classes=zeros(natomW,natomL);
classes(occ & mod(row,2)==1)=2;
classes(occ & ((mod(row,4)==2 & mod(col,2)==1) | (mod(row,4)==0 & mod(col,2)==0)))=3;
classes(occ & classes==0)=1;state.classes=classes;
names={'hex_constant','hdmi_constant','hdmi_constant1','hdmi_constant2','hdmi_constant3'};
for k=1:numel(names),state.lattice.(names{k})=eval(names{k});end
end

function restore(folder,oldpath)
cd(folder);path(oldpath);
end
