function plot_nonreciprocal(source,dest)
% Saved results only. Unnormalized means; no smoothing or solver execution.
% Angle traces mask t<0.05 ns and in-plane magnitudes below 1e-10.
assert(isfile(fullfile(source,'figure-summary.json')));
if ~isfolder(dest),mkdir(dest);end
sizes={'5x5','13x14','21x22'};counts=[16 130 336];colors=lines(6);
f=figure('Visible','off','Color','white','Position',[50 50 1550 850]);
tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');series=0;
ax=nexttile;hold(ax,'on');
labels={'low_j1_h1_s270','low_j1_h-1_s270','low_j-1_h-1_s270'};
for k=1:3
 r=loadrun(source,'2x2',labels{k});a=loadrun(source,'21x22',labels{k});
 plot(ax,a.time_ps/1000,angletrace(a,270),'Color',colors(k,:),'LineWidth',1.5);
 plot(ax,r.time_ps/1000,angletrace(r,270),'k--','LineWidth',.8,'HandleVisibility','off');series=series+2;
end
format(ax);ylabel(ax,'Balanced mean angle (deg)');title(ax,'Fig. 2(a–c) | N=336; dashed: reference');
legend(ax,{'+J, +100 Oe','+J, -100 Oe','-J, -100 Oe'},'Location','best');
ax=nexttile;hold(ax,'on');seeds=30:60:330;
for k=1:6
 a=loadrun(source,'21x22',sprintf('low_j1_h1_s%d',seeds(k)));
 plot(ax,a.time_ps/1000,angletrace(a,seeds(k)),'Color',colors(k,:));series=series+1;
end
format(ax);ylabel(ax,'Balanced mean angle (deg)');title(ax,'Fig. 2(d) | N=336; six starts');
legend(ax,compose('initial %g°',mod(60-seeds+180,360)-180),'Location','best','FontSize',8);
ax=nexttile;hold(ax,'on');fields=[175:25:400 450 500];
for q=1:4
 allsizes={'2x2',sizes{:}};values=zeros(size(fields));
 for k=1:numel(fields)
  a=loadrun(source,allsizes{q},sprintf('field%d',fields(k)));
  values(k)=mod(atan2d(a.balanced_y(end),a.balanced_z(end)),360);
 end
 styles={'k--o','-s','-^','-x'};plot(ax,fields,values,styles{q},'LineWidth',1.1);series=series+1;
end
grid(ax,'on');xlabel(ax,'Field (Oe)');ylabel(ax,'Final balanced angle (deg)');title(ax,'Fig. 5(a) | sampled transition');
legend(ax,{'Reference','N=16','N=130','N=336'},'Location','best');
for panel=1:2
 ax=nexttile;hold(ax,'on');prefixes={'balanced','site'};prefix=prefixes{panel};
 for q=1:3
  a=loadrun(source,sizes{q},'low_j1_h1_s210');
  plot(ax,a.time_ps/1000,a.([prefix '_z']),'Color',colors(q,:),'LineWidth',1.3);series=series+1;
 end
 r=loadrun(source,'2x2','low_j1_h1_s210');plot(ax,r.time_ps/1000,r.balanced_z,'k--','LineWidth',1);series=series+1;
 format(ax);ylabel(ax,[prefix ' mean z']);title(ax,['Fig. 4(a) | ' prefix ' mean']);
 legend(ax,{'N=16','N=130','N=336','Reference'},'Location','best');
end
ax=nexttile;hold(ax,'on');
for q=1:3
 path=fullfile(source,sprintf('%s-current-low_j1_h1_s210-dt5-p0.001-observables.csv',sizes{q}));
 a=readtable(path);validate(a);good=a.max_sublattice_spread>0;
 semilogy(ax,a.time_ps(good)/1000,a.max_sublattice_spread(good),'Color',colors(q,:));series=series+1;
end
set(ax,'YScale','log');format(ax);ylabel(ax,'Maximum within-sublattice spread');
title(ax,'Perturbation at t=0 | amplitude 0.001');legend(ax,compose('N=%d',counts),'Location','best');
sgtitle(f,'Nonreciprocal exchange/DMI with uniform anisotropy | pulses from 1 to 4 ns');
assert(series==27);drawnow;
exportgraphics(f,fullfile(dest,'nonreciprocal-switching.png'),'Resolution',150);
fout=fopen(fullfile(dest,'plot-checks.json'),'w');assert(fout>=0);
fprintf(fout,'%s\n',jsonencode(struct('series',series,'samples_per_trace',20001,'smoothing',false,'angle_start_ns',.05,'undefined_threshold',1e-10)));fclose(fout);
close(f);
end
function t=loadrun(source,size,label)
t=readtable(fullfile(source,sprintf('%s-current-%s-dt5-p0-observables.csv',size,label)));validate(t);
end
function validate(t)
assert(height(t)==20001 && max(abs(t.time_ps-(0:.25:5000)'))<1e-8);
assert(all(isfinite(t.balanced_y)) && all(isfinite(t.balanced_z)));
end
function phi=angletrace(t,seed)
idx=find(t.time_ps>=50 & hypot(t.balanced_y,t.balanced_z)>1e-10);phi=nan(height(t),1);
if isempty(idx),return;end
starts=[1;find(diff(idx)>1)+1];ends=[starts(2:end)-1;numel(idx)];
for j=1:numel(starts)
 q=starts(j):ends(j);v=unwrap(atan2(t.balanced_y(idx(q)),t.balanced_z(idx(q))))*180/pi;
 target=mod(60-seed+180,360)-180;phi(idx(q))=v+360*round((target-v(1))/360);
end
end
function format(ax)
xlim(ax,[0 5]);xlabel(ax,'Time (ns)');grid(ax,'on');set(ax,'FontSize',10,'Box','on');
xline(ax,1,':','HandleVisibility','off');xline(ax,4,':','HandleVisibility','off');
end
