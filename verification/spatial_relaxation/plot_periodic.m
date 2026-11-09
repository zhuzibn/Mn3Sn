function plot_periodic(out)
% Saved-data diagnostics only. Means are not normalized; no smoothing/decimation.
S=jsondecode(fileread(fullfile(out,'periodic-summary.json')));
checks=jsondecode(fileread(fullfile(out,'periodic-checks.json')));
phis=30:60:330;sizes={'8x8','12x12','20x20'};counts=[48 108 300];colors=lines(6);
fig=figure('Visible','off','Color','w','Position',[40 40 1500 1000]);
cleanup=onCleanup(@() close(fig));layout=tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
series=0;omitted=0;floor_count=0;
for q=1:3
 for row=1:3
  ax=nexttile(layout,(row-1)*3+q);hold(ax,'on');
  for k=1:6
   tag=sprintf('%s-phi%d-perturbed',sizes{q},phis(k));
   a=readtable(fullfile(out,[tag '-angles.csv']));assert(height(a)==1001);
   assert(max(abs(a.time_ns-(0:1000)'*.005))<1e-12);
   if row==1
    y=a.phi_deg;mask=a.time_ns<.05 | a.inplane_mean<1e-10 | ~isfinite(y);
    y(mask)=NaN;valid=~mask;y(valid)=rad2deg(unwrap(deg2rad(y(valid))));
    first=find(valid,1);assert(~isempty(first));y=y+360*round((phis(k)-y(first))/360);
    omitted=omitted+sum(mask);
   elseif row==2
    y=a.max_torque_T;floor_count=floor_count+sum(y<1e-15);y=max(y,1e-15);
   else
    y=a.sublattice_spread;floor_count=floor_count+sum(y<1e-15);y=max(y,1e-15);
   end
   h=plot(ax,a.time_ns,y,'Color',colors(k,:),'LineWidth',1.3,'DisplayName',sprintf('%g deg',phis(k)));
   assert(isequaln(h.XData(:),a.time_ns) && isequaln(h.YData(:),y));series=series+1;
   r=S(strcmp({S.label},tag));assert(numel(r)==1);
   assert(abs(mod(a.phi_deg(end)-r.final_phi_deg+180,360)-180)<1e-10);
  end
  xlim(ax,[0 5]);grid(ax,'on');set(ax,'FontSize',11);xlabel(ax,'Time (ns)');
  if row==1
   ylim(ax,[0 360]);yticks(ax,phis);ylabel(ax,'Mean-spin angle (deg)');
   title(ax,sprintf('%s periodic array / %d spins',sizes{q},counts(q)));
   if q==3,legend(ax,'Location','eastoutside');end
  elseif row==2
   set(ax,'YScale','log');ylim(ax,[1e-14 10]);yticks(ax,10.^[-14 -10 -7 -3 1]);ylabel(ax,'max |m x B_{eff}| (T)');
   yline(ax,1e-7,'k:','Stationarity threshold','HandleVisibility','off');
  else
   set(ax,'YScale','log');ylim(ax,[1e-15 .01]);yticks(ax,10.^[-15 -11 -7 -3]);ylabel(ax,'Max sublattice spread');
  end
 end
end
title(layout,{'Periodic relaxation: six reference seeds with spatial perturbation 0.001', ...
 'Reciprocal J/2, D/2 · uniform K · no drives · 5 fs steps · 5 ps saved samples'},'FontSize',16);
exportgraphics(fig,fullfile(out,'periodic-relaxation.png'),'Resolution',150);
clear cleanup;

fig=figure('Visible','off','Color','w','Position',[40 40 1250 760]);
cleanup=onCleanup(@() close(fig));layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
spincolors=[0 .447 .698;.835 .369 0;0 .62 .45];base=[.55 -.3;0 .65;-.55 -.3];
for k=1:6
 r=S(strcmp({S.tag},'20x20') & strcmp({S.kind},'perturbed') & [S.phi]==phis(k));assert(numel(r)==1);
 m=r.sublattice_spins;ax=nexttile(layout);hold(ax,'on');axis(ax,'equal');
 xlim(ax,[-1.5 1.5]);ylim(ax,[-2.15 1.6]);plot(ax,base([1 2 3 1],1),base([1 2 3 1],2),'-','Color',[.8 .8 .8]);
 for j=1:3
  plot(ax,base(j,1),base(j,2),'o','MarkerFaceColor',spincolors(j,:),'Color',spincolors(j,:));
  quiver(ax,base(j,1),base(j,2),.7*m(j,2),.7*m(j,3),0,'Color',spincolors(j,:),'LineWidth',2,'MaxHeadSize',.45);
 end
 net=mean(m,1);quiver(ax,0,-1.05,180*net(2),180*net(3),0,'k','LineWidth',1.6,'MaxHeadSize',.4);
 title(ax,sprintf('%g deg  ->  %.6f deg',phis(k),r.final_phi_deg),'FontSize',13);
 text(ax,0,-1.9,sprintf('E/N = %.9f meV',r.energy_meV_per_spin),'HorizontalAlignment','center');
 set(ax,'XTick',[],'YTick',[],'XColor','none','YColor','none');
end
title(layout,{'Six relaxed states on the 300-spin periodic array', ...
 'a blue / b orange / c green · black mean spin x180 · +z up, +y right · schematic positions'},'FontSize',15);
exportgraphics(fig,fullfile(out,'periodic-six-states.png'),'Resolution',150);
plotchecks=struct('trajectory_series',series,'samples_per_series',1001, ...
 'excluded_mean_angle_samples',omitted,'log_values_floored_at_1e_15',floor_count, ...
 'mean_normalization','none','smoothing','none','simulation_rerun',false);
f=fopen(fullfile(out,'plot-checks.json'),'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(plotchecks,PrettyPrint=true));fclose(f);
save(fullfile(out,'periodic-results.mat'),'S','checks','plotchecks');
fprintf('Rendered periodic trajectories and final states; %d checked line series.\n',series);
end
