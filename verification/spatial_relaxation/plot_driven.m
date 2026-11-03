function plot_driven(source,dest)
% Saved candidate/reference trajectories only; no simulations or smoothing.
% Angles omit t<0.05 ns and undefined in-plane means <1e-10.
% Mean components retain startup and are never normalized.
assert(isfile(fullfile(source,'figure-summary.json')));
sizes={'5x5','13x14','21x22'};counts=[16 130 336];seeds=30:60:330;
colors=lines(6);figures={};series=0;
f=figure('Visible','off','Color','white','Position',[50 50 1600 1050]);
tiledlayout(f,3,4,'TileSpacing','compact','Padding','compact');
labels={'low_j1_h1_s270','low_j1_h-1_s270','low_j-1_h-1_s270'};
titles={'(a) +J, +100 Oe','(b) +J, -100 Oe','(c) -J, -100 Oe','(d) +J, +100 Oe; six starts'};
for row=1:3
 for col=1:4
  ax=nexttile;hold(ax,'on');
  if col<4
   ref=loadrun(source,'2x2','current',labels{col});candidate=loadrun(source,sizes{row},'triangle',labels{col});
   plot(ax,ref.time_ps/1000,traceangle(ref,270),'--','Color',[.25 .25 .25],'LineWidth',1.1);
   plot(ax,candidate.time_ps/1000,traceangle(candidate,270),'-','Color',[0 .45 .7],'LineWidth',1.3);
   series=series+2;
   if row==1,legend(ax,{'Three-spin reference','Unequal-anisotropy array'},'Location','best','FontSize',8);end
  else
   handles=gobjects(1,6);
   for k=1:6
    label=sprintf('low_j1_h1_s%d',seeds(k));ref=loadrun(source,'2x2','current',label);candidate=loadrun(source,sizes{row},'triangle',label);
    plot(ax,ref.time_ps/1000,traceangle(ref,seeds(k)),'--','Color',[.7 .7 .7],'LineWidth',.7,'HandleVisibility','off');
    handles(k)=plot(ax,candidate.time_ps/1000,traceangle(candidate,seeds(k)),'Color',colors(k,:),'LineWidth',1);
    series=series+2;
   end
   if row==1,legend(ax,handles,compose('initial %g°',mod(60-seeds+180,360)-180),'Location','best','FontSize',8);end
  end
  format(ax);title(ax,sprintf('%s | N=%d',titles{col},counts(row)),'FontSize',10);
  ylabel(ax,'Balanced mean angle (deg)');
 end
end
sgtitle(f,'Figure 2 comparison | identical drive amplitudes and 1–4 ns pulses');
exportgraphics(f,fullfile(dest,'anisotropy-switching-fig2.png'),'Resolution',150);figures{end+1}=f;

f=figure('Visible','off','Color','white','Position',[50 50 1450 800]);
tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
ref=loadrun(source,'2x2','current','low_j1_h1_s210');cc=[.85 .15 .15;.1 .6 .3;0 .35 .85];
for row=1:2
 for col=1:3
  ax=nexttile;hold(ax,'on');candidate=loadrun(source,sizes{col},'triangle','low_j1_h1_s210');
  prefixes={'balanced','site'};prefix=prefixes{row};handles=gobjects(1,3);
  for k=1:3
   axis='xyz';field=[prefix '_' axis(k)];
   if row==1,plot(ax,ref.time_ps/1000,ref.(['balanced_' axis(k)]),'--','Color',cc(k,:),'LineWidth',1,'HandleVisibility','off');series=series+1;end
   handles(k)=plot(ax,candidate.time_ps/1000,candidate.(field),'Color',cc(k,:),'LineWidth',1.2);series=series+1;
  end
  format(ax);yline(ax,0,':','HandleVisibility','off');ylabel(ax,'Unnormalized mean spin');
  title(ax,sprintf('%s mean | N=%d',prefix,counts(col)));legend(ax,handles,{'m_x','m_y','m_z'},'Location','best');
 end
end
sgtitle(f,'Figure 4(a) comparison | +J, +100 Oe | dashed: three-spin reference');
exportgraphics(f,fullfile(dest,'anisotropy-switching-fig4a.png'),'Resolution',150);figures{end+1}=f;

f=figure('Visible','off','Color','white','Position',[50 50 1450 1100]);
tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact');
fields=[175:25:400 450 500];
for row=1:3
 for col=1:3
  ax=nexttile;hold(ax,'on');
  if row==1
   actual=zeros(size(fields));reference=actual;
   for k=1:numel(fields)
    label=sprintf('field%d',fields(k));candidate=loadrun(source,sizes{col},'triangle',label);ref=loadrun(source,'2x2','current',label);
    actual(k)=mod(atan2d(candidate.balanced_y(end),candidate.balanced_z(end)),360);
    reference(k)=mod(atan2d(ref.balanced_y(end),ref.balanced_z(end)),360);
   end
   plot(ax,fields,reference,'k--o','LineWidth',1.1);plot(ax,fields,actual,'-s','Color',[0 .45 .7],'LineWidth',1.3);series=series+2;
   xlabel(ax,'Field (Oe)');ylabel(ax,'Final balanced angle (deg)');title(ax,sprintf('(a) Field sweep | N=%d',counts(col)));
   limits=[floor(min([actual reference])/30)*30-15 ceil(max([actual reference])/30)*30+15];
   ylim(ax,limits);grid(ax,'on');legend(ax,{'Reference','Candidate'},'Location','best');
  else
   sign='p';if row==3,sign='m';end
   handles=gobjects(1,6);
   for k=1:6
    label=sprintf('high%s%d',sign,seeds(k));candidate=loadrun(source,sizes{col},'triangle',label);ref=loadrun(source,'2x2','current',label);
    plot(ax,ref.time_ps/1000,traceangle(ref,seeds(k)),'--','Color',[.7 .7 .7],'LineWidth',.7,'HandleVisibility','off');
    handles(k)=plot(ax,candidate.time_ps/1000,traceangle(candidate,seeds(k)),'Color',colors(k,:),'LineWidth',1);series=series+2;
   end
   format(ax);ylabel(ax,'Balanced mean angle (deg)');
   current='+J';letter='b';if row==3,current='-J';letter='c';end
   title(ax,sprintf('(%s) %s, +1000 Oe | N=%d',letter,current,counts(col)));
   if col==3,legend(ax,handles,compose('initial %g°',mod(60-seeds+180,360)-180),'Location','best','FontSize',8);end
  end
 end
end
sgtitle(f,'Figure 5 comparison | gray dashed: full-strength three-spin reference');
exportgraphics(f,fullfile(dest,'anisotropy-switching-fig5.png'),'Resolution',150);figures{end+1}=f;

f=figure('Visible','off','Color','white','Position',[50 50 1200 450]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for k=1:2
 labels={'low_j1_h1_s210','field275'};ax=nexttile;hold(ax,'on');
 ref=loadrun(source,'2x2','current',labels{k});candidate=loadrun(source,'2x2','triangle',labels{k});
 plot(ax,ref.time_ps/1000,traceangle(ref,210),'k--','LineWidth',1.2);
 plot(ax,candidate.time_ps/1000,traceangle(candidate,210),'Color',[.8 .3 .1],'LineWidth',1.2);series=series+2;
 format(ax);ylabel(ax,'Mean angle (deg)');fields2=[100 275];title(ax,sprintf('+J, +%d Oe',fields2(k)));
 legend(ax,{'Full-strength reference','Half-strength candidate triangle'},'Location','best');
end
sgtitle(f,'Three-spin control | candidate J/2, D/2, K/2 with unchanged external drives');
exportgraphics(f,fullfile(dest,'anisotropy-switching-triangle.png'),'Resolution',150);figures{end+1}=f;
assert(series==163);
fout=fopen(fullfile(source,'plot-checks.json'),'w');assert(fout>=0);
fprintf(fout,'%s\n',jsonencode(struct('figures',4,'series',series,'samples_per_trace',20001,'smoothing',false,'angle_start_ns',.05,'undefined_threshold',1e-10)));fclose(fout);
for k=1:numel(figures),close(figures{k});end
end
function t=loadrun(source,size,mode,label)
path=fullfile(source,sprintf('%s-%s-%s-dt5-p0-observables.csv',size,mode,label));
t=readtable(path);assert(height(t)==20001 && max(abs(t.time_ps-(0:.25:5000)'))<1e-8);
assert(all(isfinite(t.balanced_y)) && all(isfinite(t.balanced_z)));
end
function phi=traceangle(t,seed)
good=t.time_ps>=50 & hypot(t.balanced_y,t.balanced_z)>1e-10;
phi=nan(height(t),1);idx=find(good);raw=atan2(t.balanced_y(idx),t.balanced_z(idx));
if isempty(idx),return;end
% Unwrap each contiguous valid interval, never bridging an undefined mean.
starts=[1;find(diff(idx)>1)+1];ends=[starts(2:end)-1;numel(idx)];
for j=1:numel(starts)
 q=starts(j):ends(j);v=unwrap(raw(q))*180/pi;
 target=mod(60-seed+180,360)-180;v=v+360*round((target-v(1))/360);phi(idx(q))=v;
end
end
function format(ax)
xlim(ax,[0 5]);xlabel(ax,'Time (ns)');grid(ax,'on');set(ax,'FontSize',10,'Box','on');
xline(ax,1,':','HandleVisibility','off');xline(ax,4,':','HandleVisibility','off');
end
