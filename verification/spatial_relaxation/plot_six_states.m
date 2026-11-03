function plot_six_states(out,destination,phase)
% Saved trajectories only; angles from +z toward +y, modulo 360 degrees.
S=jsondecode(fileread(fullfile(out,'six-state-summary.json')));
phis=30:60:330;
if strcmp(phase,'reference')
 fig=figure('Visible','off','Color','w','Position',[40 40 1200 750]);
 cleanup=onCleanup(@() close(fig));layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
 for k=1:6
  r=S(strcmp({S.size},'2x2') & [S.initial_phi_deg]==phis(k) & [S.perturb]==0);assert(numel(r)==1);
  assert(abs(mod(r.final_phi_deg-phis(k)+180,360)-180)<.1);
  ax=nexttile(layout);drawstate(ax,r.sublattice_spins);
  title(ax,sprintf('State %c: phi = %g degrees','a'+k-1,phis(k)));
  text(ax,0,-1.85,sprintf('spin angles: %.2f, %.2f, %.2f deg',mod(atan2d(r.sublattice_spins(:,2),r.sublattice_spins(:,3)),360)),...
   'HorizontalAlignment','center','FontSize',10);
 end
 title(layout,{'Six relaxed three-sublattice reference states (zero field and current)',...
  'Blue: spin a   Orange: spin b   Green: spin c   Black: mean spin x 180 | +z up, +y right'},'FontSize',15);
 exportgraphics(fig,fullfile(destination,'six-reference-states.png'),'Resolution',150);
 fprintf('Rendered six reference states from actual relaxed vectors.\n');return;
end
assert(strcmp(phase,'comparison'));
fig=figure('Visible','off','Color','w','Position',[40 40 1400 1150]);
cleanup=onCleanup(@() close(fig));layout=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
colors=lines(6);modes={'current','reciprocal'};labels={'Nonreciprocal weights','Reciprocal J/2, D/2 bonds'};
for q=1:2
 ax=nexttile(layout,q);hold(ax,'on');
 for k=1:6
  a=readtable(fullfile(out,sprintf('21x22-%s-phi%d-p0.001-angles.csv',modes{q},phis(k))));
  angle=rad2deg(unwrap(deg2rad(a.order_phi_deg)));
  angle=angle+360*round((phis(k)-angle(1))/360);
  plot(ax,a.time_ns*1000,angle,'Color',colors(k,:),'LineWidth',1.35,'DisplayName',sprintf('%c: %g deg','a'+k-1,phis(k)));
 end
 for phi=phis,yline(ax,phi,':','Color',[.7 .7 .7],'HandleVisibility','off');end
 xlabel(ax,'Time (ps)');ylabel(ax,'Spin-a orientation: 60 - theta_a (deg)');
 xlim(ax,[0 150]);ylim(ax,[-30 390]);yticks(ax,phis);grid(ax,'on');title(ax,labels{q});set(ax,'FontSize',11);
 legend(ax,'Location','eastoutside');
 ax=nexttile(layout,q+2);hold(ax,'on');
 for k=1:6
  a=readtable(fullfile(out,sprintf('21x22-%s-phi%d-p0.001-angles.csv',modes{q},phis(k))));
  % Undefined/compensated and initial transient moments are excluded explicitly.
  angle=a.balanced_phi_deg;angle(a.balanced_magnitude<1e-10 | a.time_ns<.05)=NaN;
  valid=isfinite(angle);angle(valid)=rad2deg(unwrap(deg2rad(angle(valid))));
  first=find(valid,1);angle=angle+360*round((phis(k)-angle(first))/360);
  plot(ax,a.time_ns,angle,'Color',colors(k,:),'LineWidth',1.35,'DisplayName',sprintf('%c: %g deg','a'+k-1,phis(k)));
 end
 for phi=phis,yline(ax,phi,':','Color',[.7 .7 .7],'HandleVisibility','off');end
 xlabel(ax,'Time (ns)');ylabel(ax,'Balanced three-sublattice mean angle (deg)');
 xlim(ax,[0 5]);ylim(ax,[-30 390]);yticks(ax,phis);grid(ax,'on');title(ax,'Mean-angle destinations over the full run');set(ax,'FontSize',11);
 legend(ax,'Location','eastoutside');
 ax=nexttile(layout,q+4);hold(ax,'on');
 for k=1:6
  a=readtable(fullfile(out,sprintf('21x22-%s-phi%d-p0.001-angles.csv',modes{q},phis(k))));
  semilogy(ax,a.time_ns,a.max_torque_T,'Color',colors(k,:),'LineWidth',1.2);
 end
 set(ax,'YScale','log','FontSize',11);xlim(ax,[0 5]);ylim(ax,[1e-14 100]);grid(ax,'on');
 yline(ax,1e-6,'k:','Stationarity criterion');xlabel(ax,'Time (ns)');ylabel(ax,'Maximum |m x B_{eff}| (T)');
 title(ax,'Residual torque for the same six trajectories');
end
title(layout,{'Six-state test on the 21 x 22 open array (336 spins)',...
 'Six reference seeds + spatial perturbation 0.001 | alpha = 0.003 | dt = 5 fs | no drives'},'FontSize',16);
exportgraphics(fig,fullfile(destination,'six-state-relaxation-comparison.png'),'Resolution',150);
clear cleanup;

fig=figure('Visible','off','Color','w','Position',[40 40 1600 650]);
cleanup=onCleanup(@() close(fig));layout=tiledlayout(fig,2,6,'TileSpacing','compact','Padding','compact');
for q=1:2
 for k=1:6
  r=S(strcmp({S.size},'21x22') & strcmp({S.mode},modes{q}) & [S.initial_phi_deg]==phis(k) & [S.perturb]==.001);assert(numel(r)==1);
  ax=nexttile(layout);drawstate(ax,r.sublattice_spins);
  title(ax,sprintf('%c: %g -> %.2f deg','a'+k-1,phis(k),r.final_phi_deg),'FontSize',11);
  if k==1,ylabel(ax,labels{q},'FontSize',12);end
 end
end
title(layout,{'Relaxed sublattice-average spins for each of the six starts',...
 'Top: nonreciprocal | Bottom: reciprocal | Blue a, orange b, green c | black balanced mean x 180'},'FontSize',15);
exportgraphics(fig,fullfile(destination,'six-state-final-spins.png'),'Resolution',150);
fprintf('Rendered comparison trajectories and twelve final spin diagrams.\n');
end

function drawstate(ax,m)
hold(ax,'on');axis(ax,'equal');xlim(ax,[-1.5 1.5]);ylim(ax,[-2 1.6]);
base=[.55 -.3;0 .65;-.55 -.3];colors=[0 .447 .698;.835 .369 0;0 .62 .45];
plot(ax,base([1 2 3 1],1),base([1 2 3 1],2),'-','Color',[.8 .8 .8]);
for j=1:3
 plot(ax,base(j,1),base(j,2),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',5);
 quiver(ax,base(j,1),base(j,2),.7*m(j,2),.7*m(j,3),0,'Color',colors(j,:),'LineWidth',2,'MaxHeadSize',.45);
end
meanm=mean(m,1);quiver(ax,0,-1.05,180*meanm(2),180*meanm(3),0,'k','LineWidth',1.6,'MaxHeadSize',.4);
set(ax,'XTick',[],'YTick',[],'Box','off','XColor','none','YColor','none');
end
