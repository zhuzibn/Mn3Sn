function p=mn3sn_protocol(name,value)
% One relaxation or PRB 109, 134433 (2024) figure trajectory.
% For fig2d/fig5b/fig5c, value is the first-spin seed in degrees.
% For fig5a, value is the field in Oe. Other presets take no value.
if nargin<1, name='relaxation'; end
if nargin<2, value=[]; end
name=validatestring(name,{'relaxation','fig2a','fig2b','fig2c','fig2d','fig4a','fig5a','fig5b','fig5c'});
p=struct('name',name,'natomW',2,'natomL',2,'seed_deg',270, ...
    'current',5e10,'field',[0 .01 0],'pulse_start',1e-9,'pulse_end',4e-9, ...
    'runtime',5e-9,'tstep',5e-15,'gpusave',10e-12,'savetstep',50);
if ismember(name,{'fig2d','fig5b','fig5c'})
    if isempty(value),value=210;end
    validateattributes(value,{'numeric'},{'scalar','real','finite'});
    assert(ismember(value,30:60:330),'Mn3Sn:Seed','Choose a seed from 30:60:330 degrees.');
    p.seed_deg=value;
elseif strcmp(name,'fig5a')
    if isempty(value),value=250;end
    validateattributes(value,{'numeric'},{'scalar','real','finite'});
    assert(ismember(value,[175 200 225 250 275 300 325 350 375 400 450 500]), ...
        'Mn3Sn:Field','Choose a tested Fig. 5(a) field: 175:25:400, 450 or 500 Oe.');
else
    assert(isempty(value),'Mn3Sn:ProtocolValue','This preset takes no second argument.');
end
switch name
    case 'relaxation'
        p.natomW=5;p.natomL=5;p.seed_deg=121;
        p.current=0;p.field=[0 0 0];p.pulse_start=0;p.pulse_end=Inf;
        p.runtime=10e-12;p.tstep=.5e-15;p.gpusave=1e-12;p.savetstep=200;
    case 'fig2b'
        p.field=[0 -.01 0];
    case 'fig2c'
        p.current=-5e10;p.field=[0 -.01 0];
    case 'fig4a'
        p.seed_deg=210;
    case 'fig5a'
        p.seed_deg=210;p.field=[0 value*1e-4 0];
    case 'fig5b'
        p.field=[0 .1 0];
    case 'fig5c'
        p.current=-5e10;p.field=[0 .1 0];
end
end
