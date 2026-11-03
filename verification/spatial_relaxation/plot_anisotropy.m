function plot_anisotropy(previous,out,destination)
% Saved 336-spin data only; no simulation or modifications to raw inputs.
% Balanced mean angles omit t<0.05 ns / |mean_yz|<1e-10, unwrap from each seed.
% Final arrows use actual sublattice averages, black balanced mean amplified 180x.
S=jsondecode(fileread(fullfile(out,'six-state-summary.json')));
checks=jsondecode(fileread(fullfile(out,'trial-checks.json')));assert(checks.main_trajectories==60);
phis=30:60:330;colors=lines(6);styles={'-','--',':','-','--',':'};
fig=figure('Visible','off','Color','w','Position',[40 40 1400 800]);
cleanup=onCleanup(@() close(fig));layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
folders={previous,out};modes={'reciprocal','triangle'};
labels={'Uniform K: two destinations','K_i = K z_i / 4: six destinations'};
for q=1:2
 ax=nexttile(layout,q);hold(ax,'on');
 torqueax=nexttile(layout,q+2);hold(torqueax,'on');
 for k=1:6
  T=readtable(fullfile(folders{q},sprintf('21x22-%s-phi%d-p0.001-angles.csv',modes{q},phis(k))));
  assert(height(T)==1001 && abs(T.time_ns(end)-5)<1e-12);
  v=T.balanced_phi_deg;v(T.time_ns<.05 | T.balanced_magnitude<1e-10)=NaN;
  valid=isfinite(v);v(valid)=rad2deg(unwrap(deg2rad(v(valid))));
  first=find(valid,1);v=v+360*round((phis(k)-v(first))/360);
  assert(max(abs(mod(v(valid)-T.balanced_phi_deg(valid)+180,360)-180))<1e-10);
  plot(ax,T.time_ns,v,styles{k},'Color',colors(k,:),'LineWidth',1.5,'DisplayName',sprintf('%g deg',phis(k)));
  semilogy(torqueax,T.time_ns,T.max_torque_T,styles{k},'Color',colors(k,:),'LineWidth',1.2);
 end
 assert(numel(findobj(ax,'Type','line'))==6 && numel(findobj(torqueax,'Type','line'))==6);
 for phi=phis,yline(ax,phi,':','Color',[.8 .8 .8],'HandleVisibility','off');end
 title(ax,labels{q},'FontSize',15);xlabel(ax,'Time (ns)');ylabel(ax,'Balanced mean angle (deg)');
 set(ax,'FontSize',11);xlim(ax,[0 5]);ylim(ax,[0 360]);yticks(ax,phis);grid(ax,'on');
 legend(ax,'Location','eastoutside');
 set(torqueax,'YScale','log','FontSize',11);xlim(torqueax,[0 5]);ylim(torqueax,[1e-14 10]);grid(torqueax,'on');
 yticks(torqueax,10.^(-12:3:0));
 yline(torqueax,1e-6,'k:','Stationarity threshold');xlabel(torqueax,'Time (ns)');ylabel(torqueax,'Maximum |m x B| (T)');
 title(torqueax,'Residual torque for the same six starts','FontSize',13);
end
title(layout,{'Reciprocal J/2, D/2 on the same 336-spin open array',...
 'Six reference seeds + 0.001 spatial perturbation | dt = 5 fs | no drives'},'FontSize',16);
drawnow;exportgraphics(fig,fullfile(destination,'anisotropy-six-state.png'),'Resolution',150);
clear cleanup;
fig=figure('Visible','off','Color','w','Position',[40 40 1200 750]);
cleanup=onCleanup(@() close(fig));layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for k=1:6
 r=S(strcmp({S.size},'21x22') & [S.initial_phi_deg]==phis(k) & [S.perturb]==.001);assert(numel(r)==1);
 assert(abs(mod(r.final_phi_deg-phis(k)+180,360)-180)<.1);
 ax=nexttile(layout,k);drawstate(ax,r.sublattice_spins);
 title(ax,sprintf('%g -> %.2f deg',phis(k),r.final_phi_deg),'FontSize',12);
end
title(layout,{'Six final states with reciprocal J/2, D/2 and K_i = K z_i / 4',...
 '336 spins | a blue, b orange, c green | balanced mean x 180 black | +z up, +y right'},'FontSize',15);
drawnow;exportgraphics(fig,fullfile(destination,'anisotropy-final-spins.png'),'Resolution',150);
fprintf('Exported anisotropy comparison: 24 traces and six final spin configurations.\n');
end

function drawstate(ax,m)
hold(ax,'on');axis(ax,'equal');xlim(ax,[-1.5 1.5]);ylim(ax,[-1.9 1.6]);
base=[.55 -.3;0 .65;-.55 -.3];colors=[0 .447 .698;.835 .369 0;0 .62 .45];
plot(ax,base([1 2 3 1],1),base([1 2 3 1],2),'-','Color',[.8 .8 .8]);
for j=1:3
 plot(ax,base(j,1),base(j,2),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',5);
 quiver(ax,base(j,1),base(j,2),.7*m(j,2),.7*m(j,3),0,'Color',colors(j,:),'LineWidth',2,'MaxHeadSize',.45);
end
meanm=mean(m,1);quiver(ax,0,-1.05,180*meanm(2),180*meanm(3),0,'k','LineWidth',1.6,'MaxHeadSize',.4);
set(ax,'XTick',[],'YTick',[],'Box','off','XColor','none','YColor','none');
end
