% This script draws the histogram for rho values (Figure 1.(a))
% It has to be run after MAIN.m, while the variable rho is still in the
% workspace or by loading it from saved data

ftsz = 25;

histogram(rho,'FaceColor','b','FaceAlpha',0.7,'EdgeAlpha',1,'LineWidth',2)
hold on
grid on
intr = 'latex';

ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = ftsz;

xlabel('$\hat{\rho}$','fontsize',ftsz,'Interpreter',intr)
maxrep = 40;
yticks(0:5:maxrep)
xticks(1:1:25)
ylim([0 maxrep])
xlim([0 25])

drawnow