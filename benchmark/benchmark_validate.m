function o=benchmark_validate(r)
%BENCHMARK_VALIDATE Recompute invariants/observables; reject malformed results.
assert(isstruct(r) && isscalar(r) && r.schema_version==1,'Unsupported schema.');
assert(ismember(r.mode,{'quick','current'}) && ...
    ismember(r.case_name,{'relaxation','driven_nonuniform'}),'Unknown case/mode.');
c=r.config;W=4+4*strcmp(r.mode,'current');
assert(c.natomW==W && c.natomL==W && c.rk4==1 && c.bc==0);
assert(isequal(size(r.atomtype),[W W]) && all(ismember(r.atomtype(:),[1 3])));
occ=r.atomtype==1;empty=~occ;N=nnz(occ);assert(N==3*W*W/4);
assert(isequal(size(r.classes),[W W]) && all(r.classes(empty)==0));
assert(all(arrayfun(@(q)nnz(r.classes==q & occ)==N/3,1:3)));
assert(isequal(size(r.axes),[W W 3]) && isreal(r.axes) && all(isfinite(r.axes),'all'));
assert(isequal(size(r.initial),[W W 3]) && isequal(size(r.generated_initial),[W W 3]));
assert(isreal(r.initial) && all(isfinite(r.initial),'all') && ...
    isreal(r.generated_initial) && all(isfinite(r.generated_initial),'all'));
if strcmp(r.case_name,'relaxation'),updates=20000;frames=101;active=20000;else,updates=60;frames=61;active=31;end
assert(r.actual_updates==updates && r.pulse_active_updates==active,'Update count mismatch.');
assert(round(c.runtime/c.tstep)==updates && ...
    c.gpurun_number*(c.gpusteps-1)==updates,'Configuration clock mismatch.');
steps=unique([0:c.savetstep:updates updates]);
assert(numel(steps)==frames && isequal(r.saved_steps,steps),'Saved step mismatch.');
assert(isequal(size(r.t),size(steps)) && isreal(r.t) && all(isfinite(r.t)));
assert(max(abs(r.t-steps*c.tstep))<=8*eps(c.runtime),'Time mismatch.');
assert(sum((0:updates-1)>=round(c.pulse_start/c.tstep) & ...
    (0:updates-1)<round(c.pulse_end/c.tstep))==active,'Pulse mismatch.');
x=r.mmx;y=r.mmy;z=r.mmz;
for a={x,y,z}
    assert(isa(a{1},'double') && isreal(a{1}) && isequal(size(a{1}),[W W frames]) && ...
        all(isfinite(a{1}),'all'),'Malformed spin trajectory.');
    e=repmat(empty,1,1,frames);assert(all(a{1}(e)==0),'Nonzero vacancy.');
end
assert(isequal(cat(3,x(:,:,1),y(:,:,1),z(:,:,1)),r.initial),'Initial frame mismatch.');
n=sqrt(x.^2+y.^2+z.^2);mask=repmat(occ,1,1,frames);
o.norm_error=max(abs(n(mask)-1));assert(o.norm_error<=1e-12,'Spin norm failure.');
o.mean=zeros(frames,3);o.class_mean=zeros(frames,3,3);
for q=1:3
    a={x,y,z};v=reshape(a{q},W*W,frames);
    o.mean(:,q)=mean(v(occ(:),:),1)';
    for cl=1:3,o.class_mean(:,q,cl)=mean(v(r.classes(:)==cl,:),1)';end
end
o.vacancy_max=0;
end
