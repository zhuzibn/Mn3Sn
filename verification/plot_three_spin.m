function plot_three_spin(out)
% Plot saved data only. No solver invocation or altered numerical parameters.
% Binary layout: physical seconds, ax ay az bx by bz cx cy cz.
% Signed in-plane mean angle; zero mean undefined. Break wrapping segments.
% Angle plots start at 0.05 ns to avoid the ideal compensated startup seed.
assert(jsondecode(fileread(fullfile(out,'analysis.json'))).pass);
T=readtable(fullfile(out,'summary.csv'),'VariableNamingRule','preserve');
fig=figure('Visible','off','Color','w','Position',[50 50 1250 1150]);
clean=onCleanup(@() close(fig));layout=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
colors=lines(12);fields=[175 200 225 250 275 300 325 350 375 400 450 500];
ax=nexttile(layout);hold(ax,'on');
for k=1:numel(fields)
    row=T(T.mode==2 & T.dt_fs==1.25 & T.field_Oe==fields(k),:);assert(height(row)==1);
    [time,angle]=angles(readbin(fullfile(out,[row.name{1} '.bin'])));
    if k<=7,linestyle='-';else,linestyle='--';end
    plot(ax,time,angle,'LineStyle',linestyle,'LineWidth',1.2,'Color',colors(k,:), ...
        'DisplayName',sprintf('%g Oe',fields(k)));
end
style(ax,'(a) Matched coupled RK4: field sweep, seed 210 degrees');
legend(ax,'Location','southoutside','NumColumns',6);
for panel=1:2
    ax=nexttile(layout);hold(ax,'on');current=5e10*(3-2*panel);
    for k=1:6
        seed=30+60*(k-1);
        row=T(T.mode==2 & T.dt_fs==1.25 & T.field_Oe==1000 & T.current==current & T.seed_deg==seed,:);
        assert(height(row)==1);
        [time,angle]=angles(readbin(fullfile(out,[row.name{1} '.bin'])));
        plot(ax,time,angle,'LineWidth',1.2,'Color',colors(k,:), ...
            'DisplayName',sprintf('Seed %g deg',seed));
    end
    style(ax,sprintf('(%c) 1000 Oe, current %+g A/m^2','a'+panel,current));
    legend(ax,'Location','southoutside','NumColumns',6);
end
title(layout,'Three-spin Figure 5 comparison | dt = 1.25 fs | drive 1–4 ns','FontSize',16);
exportgraphics(fig,fullfile(out,'fig5.png'),'Resolution',150);
clear clean;

fig=figure('Visible','off','Color','w','Position',[50 50 1450 1000]);
clean=onCleanup(@() close(fig));layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
colors=[0 .447 .698;.835 .369 0;.0 .62 .45;.6 .3 .65];
for mode=0:1
    ax=nexttile(layout);hold(ax,'on');
    for k=1:3
        dt=[5 1.25 .5];dt=dt(k);
        row=T(T.mode==mode & T.dt_fs==dt & T.field_Oe==300,:);assert(height(row)==1);
        [time,angle]=angles(readbin(fullfile(out,[row.name{1} '.bin'])));
        plot(ax,time,angle,'LineWidth',1.1,'Color',colors(k,:),'DisplayName',sprintf('Frozen, %g fs',dt));
    end
    row=T(T.mode==2 & T.dt_fs==1.25 & T.field_Oe==300,:);
    [time,angle]=angles(readbin(fullfile(out,[row.name{1} '.bin'])));
    plot(ax,time,angle,'k--','LineWidth',1.5,'DisplayName','Matched coupled');
    labels={'(a) Native parameters, frozen fields (300 Oe)','(b) Matched parameters, frozen fields (300 Oe)'};
    style(ax,labels{mode+1});legend(ax,'Location','southoutside','NumColumns',2);
end
ax=nexttile(layout);hold(ax,'on');
for mode=0:2
    rows=T(T.mode==mode & T.dt_fs==(0.5*(mode<2)+1.25*(mode==2)) & T.field_Oe<1000,:);
    [~,order]=sort(rows.field_Oe);rows=rows(order,:);
    labels={'Native frozen, 0.5 fs','Matched frozen, 0.5 fs','Matched coupled, 1.25 fs'};
    plot(ax,rows.field_Oe,rows.final_angle_deg,'o-','Color',colors(mode+1,:), ...
        'LineWidth',1.2,'DisplayName',labels{mode+1});
end
grid(ax,'on');xlabel(ax,'Applied field (Oe)');ylabel(ax,'Final mean angle (degrees)');
title(ax,'(c) Sampled 5 ns endpoints');legend(ax,'Location','southoutside');set(ax,'FontSize',11);
ax=nexttile(layout);hold(ax,'on');
for mode=0:2
    rows=T(T.mode==mode & T.field_Oe==300,:);[~,order]=sort(rows.dt_fs);rows=rows(order,:);
    labels={'Native frozen','Matched frozen','Matched coupled'};
    loglog(ax,rows.dt_fs,rows.max_mean_component_difference,'o-','Color',colors(mode+1,:), ...
        'LineWidth',1.2,'DisplayName',labels{mode+1});
end
R=jsondecode(fileread(fullfile(out,'refinement.json')));
for mode=0:1
    rows=R.comparisons([R.comparisons.mode]==mode);
    values=arrayfun(@(r)r.reference.max_mean_component_difference,rows);
    loglog(ax,[rows.dt_fs],values,'s:','Color',colors(mode+1,:),'LineWidth',1.2,'HandleVisibility','off');
end
yline(ax,1e-5,'k--','Target','HandleVisibility','off');grid(ax,'on');
set(ax,'XScale','log','YScale','log','FontSize',11);
xlabel(ax,'Timestep (fs)');ylabel(ax,'Maximum mean-component difference');
title(ax,'(d) Error versus reference at 1.25 fs, 300 Oe');legend(ax,'Location','southoutside');
title(layout,'Frozen fields distort dynamics despite nearly exact spin norms','FontSize',16);
exportgraphics(fig,fullfile(out,'comparison.png'),'Resolution',150);
% Verify all plotted rows originate from the checked source table.
assert(height(T)==336);
fprintf('Exported Figure 5 and comparison PNGs from %d checked trajectories.\n',height(T));
end

function style(ax,label)
grid(ax,'on');xlim(ax,[0 5]);ylim(ax,[-185 185]);yticks(ax,-180:60:180);
xlabel(ax,'Physical time (ns)');ylabel(ax,'Mean angle (degrees)');title(ax,label);
xline(ax,1,':','HandleVisibility','off');xline(ax,4,':','HandleVisibility','off');set(ax,'FontSize',11);
end
function [time,angle]=angles(data)
time=data(:,1)*1e9;my=mean(data(:,[3 6 9]),2);mz=mean(data(:,[4 7 10]),2);
angle=atan2d(my,mz);angle(hypot(my,mz)<1e-10 | time<.05)=NaN;
angle([false;abs(diff(angle))>180])=NaN;
end
function data=readbin(path)
f=fopen(path,'rb','ieee-le');assert(f>=0);cleanup=onCleanup(@() fclose(f));data=fread(f,[10 Inf],'double')';
end
