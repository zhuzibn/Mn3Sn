function plot_coupled_figures(out,plotout)
% Saved 5 fs matched-parameter trajectories only; no simulation or smoothing.
% Mean of three unit spins, signed atan2(y,z), physical seconds -> ns.
% Angle plots omit t<0.05 ns and zero means <1e-10; unwrap each valid segment.
% Optional plotout separates regenerated figures from the saved run directory.
if nargin<2,plotout=out;end
assert(jsondecode(fileread(fullfile(out,'figures-analysis.json'))).pass);
colors=lines(6);seeds=30:60:330;
fig=figure('Visible','off','Color','w','Position',[30 30 1200 850]);
cleanup=onCleanup(@() close(fig));lay=tiledlayout(fig,2,2,'TileSpacing','compact');
signs=[1 1;1 -1;-1 -1];
for panel=1:3
    ax=nexttile(lay);name=sprintf('low_j%d_h%d_s270',signs(panel,:));
    [t,~,angle]=trace(out,name);plot(ax,t,angle,'LineWidth',1.3);
    style(ax,sprintf('(%c) J %+.0e A/m^2, H %+d Oe','a'+panel-1,signs(panel,1)*5e10,signs(panel,2)*100));
end
ax=nexttile(lay);hold(ax,'on');
for k=1:6
    [t,~,angle]=trace(out,sprintf('low_j1_h1_s%d',seeds(k)));
    plot(ax,t,angle,'LineWidth',1.2,'Color',colors(k,:),'DisplayName',sprintf('Initial %g deg',wrap(60-seeds(k))));
end
style(ax,'(d) All six initial states, +J and +100 Oe');legend(ax,'Location','southoutside','NumColumns',3);
title(lay,'Figure 2 | Matched three-spin parameters | Coupled RK4, 5 fs');
angle_limits(fig);
exportgraphics(fig,fullfile(plotout,'fig2.png'),'Resolution',150);clear cleanup;

fig=figure('Visible','off','Color','w','Position',[30 30 1100 550]);
cleanup=onCleanup(@() close(fig));ax=axes(fig);hold(ax,'on');
[t,m,~]=trace(out,'low_j1_h1_s210');labels={'m_{oct,x}','m_{oct,y}','m_{oct,z}'};styles={'-','--','-.'};
for k=1:3,plot(ax,t,m(:,k),'LineWidth',1.3,'LineStyle',styles{k},'DisplayName',labels{k});end
grid(ax,'on');xlim(ax,[0 5]);xlabel(ax,'Physical time (ns)');ylabel(ax,'Mean of three unit spins');
xline(ax,1,':','HandleVisibility','off');xline(ax,4,':','HandleVisibility','off');legend(ax,'Location','best');
title(ax,{'Figure 4(a) | Mn_3Sn | +J, +100 Oe | Coupled RK4, 5 fs', ...
    'Ideal 120-degree seed: startup relaxation retained'});
exportgraphics(fig,fullfile(plotout,'fig4a.png'),'Resolution',150);clear cleanup;

fig=figure('Visible','off','Color','w','Position',[30 30 1200 1100]);
cleanup=onCleanup(@() close(fig));lay=tiledlayout(fig,3,1,'TileSpacing','compact');
ax=nexttile(lay);hold(ax,'on');fields=[175 200 225 250 275 300 325 350 375 400 450 500];cols=lines(12);
for k=1:numel(fields)
    [t,~,angle]=trace(out,sprintf('field%d',fields(k)));ls='-';if k>6,ls='--';end
    plot(ax,t,angle,'LineWidth',1.2,'LineStyle',ls,'Color',cols(k,:),'DisplayName',sprintf('%d Oe',fields(k)));
end
style(ax,'(a) Field sweep, +J, initial -150 degrees');legend(ax,'Location','southoutside','NumColumns',6);
for panel=1:2
    ax=nexttile(lay);hold(ax,'on');sgn={'p','m'};current=[5e10 -5e10];
    for k=1:6
        [t,~,angle]=trace(out,sprintf('high%s%d',sgn{panel},seeds(k)));
        plot(ax,t,angle,'LineWidth',1.2,'Color',colors(k,:),'DisplayName',sprintf('Initial %g deg',wrap(60-seeds(k))));
    end
    style(ax,sprintf('(%c) 1000 Oe, J %+.0e A/m^2','a'+panel,current(panel)));
    legend(ax,'Location','southoutside','NumColumns',6);
end
title(lay,'Figure 5 | Matched three-spin parameters | Coupled RK4, 5 fs');
angle_limits(fig);
exportgraphics(fig,fullfile(plotout,'fig5.png'),'Resolution',150);
fprintf('Exported Figures 2, 4(a), 5 from saved 5 fs runs.\n');
end

function [t,m,angle]=trace(out,label)
f=fopen(fullfile(out,[label '_m2_dt5.bin']),'rb','ieee-le');assert(f>=0);c=onCleanup(@() fclose(f));
a=fread(f,[10 Inf],'double')';assert(isequal(size(a),[20001 10]));
t=a(:,1)*1e9;m=[mean(a(:,[2 5 8]),2),mean(a(:,[3 6 9]),2),mean(a(:,[4 7 10]),2)];
angle=atan2d(m(:,2),m(:,3));angle(vecnorm(m,2,2)<1e-10 | t<.05)=NaN;
% Preserve real gaps; remove only the coordinate jump between -180 and +180.
edges=diff([false;isfinite(angle);false]);starts=find(edges==1);stops=find(edges==-1)-1;
for k=1:numel(starts)
    idx=starts(k):stops(k);angle(idx)=rad2deg(unwrap(deg2rad(angle(idx))));
end
end
function a=wrap(a),a=mod(a+180,360)-180;end
function style(ax,label)
grid(ax,'on');xlim(ax,[0 5]);
xlabel(ax,'Physical time (ns)');ylabel(ax,'Unwrapped octupole angle (deg)');title(ax,label);
xline(ax,1,':','HandleVisibility','off');xline(ax,4,':','HandleVisibility','off');
end
function angle_limits(fig)
% Common limits across comparable panels, enclosing every finite sample.
data=get(findobj(fig,'Type','line'),'YData');values=[data{:}];values=values(isfinite(values));
limits=90*[floor(min(values)/90),ceil(max(values)/90)];
for ax=reshape(findobj(fig,'Type','axes'),1,[])
    ylim(ax,limits+[-10 10]);yticks(ax,limits(1):90:limits(2));
end
end
