function plot_spatial(out,png)
% Saved scalar traces only. Time ps -> ns; no smoothing or fitted curves.
% All six current-weight curves overlap; label the representative largest one.
fig=figure('Visible','off','Color','w','Position',[100 100 1200 450]);
cleanup=onCleanup(@() close(fig));
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout);hold(ax,'on');colors=lines(5);sizes={'5x5','5x6','9x10','13x14','21x22'};
a=readtable(fullfile(out,'21x22-current-5fs-5ns.csv'));
assert(height(a)==51 && all(a.max_torque_T>0));
semilogy(ax,a.time_ps/1000,a.max_torque_T,'k--','LineWidth',1.8,'DisplayName','Current weights (336 spins)');
for k=1:5
 a=readtable(fullfile(out,[sizes{k} '-reciprocal-5fs-5ns.csv']));
 assert(height(a)==51 && all(a.max_torque_T>0));
 semilogy(ax,a.time_ps/1000,a.max_torque_T,'-o','Color',colors(k,:),'MarkerSize',3,...
  'LineWidth',1.2,'DisplayName',[strrep(sizes{k},'x',' x ') ' reciprocal']);
end
yline(ax,1e-6,':','Convergence criterion','HandleVisibility','off');
set(ax,'YScale','log','FontSize',11);ylim(ax,[1e-14 10]);xlim(ax,[0 1.5]);grid(ax,'on');
xlabel(ax,'Time (ns)');ylabel(ax,'Maximum |m x B_{eff}| (T)');title(ax,'Undriven relaxation; all runs continue to 5 ns');
legend(ax,'Location','northeast','FontSize',9);
ax=nexttile(layout);hold(ax,'on');modes={'current','reciprocal'};styles={'--','-'};
for k=1:2
 a=readtable(fullfile(out,['21x22-' modes{k} '-recovery.csv']));
 assert(height(a)==21 && all(a.max_torque_T>0));
 semilogy(ax,a.time_ps/1000,a.max_torque_T,styles{k},'Color',colors(1+3*(k-1),:),...
  'LineWidth',1.8,'DisplayName',modes{k});
end
yline(ax,1e-6,':','Convergence criterion','HandleVisibility','off');
set(ax,'YScale','log','FontSize',11);ylim(ax,[1e-14 10]);xlim(ax,[0 2]);grid(ax,'on');
xlabel(ax,'Time after perturbation (ns)');ylabel(ax,'Maximum |m x B_{eff}| (T)');
title(ax,'336 spins: recovery after ~10^{-3} rad perturbation');legend(ax,'Location','northeast');
title(layout,'Open boundaries | J, D, K matched to the three-spin reference | alpha = 0.003 | dt = 5 fs');
drawnow;exportgraphics(fig,png,'Resolution',150);
fprintf('Exported %s; 8 plotted traces checked against saved tables.\n',png);
end
