% This is a UTILITY function to draw the Figure 1.(a) in [9].
% Given
% - graphic options gra
% - the figure object fig
% - the values of beta (betas) to show on the x axis
% - the values of beta (bbar) minimizing the cost Jmin
% - the cost Jave computed by averaging the Monte Carlo experiment run
%   inside innerMC_grids
% - Jmin = Jave(bbar)
% - the x label xlab
% - the y label ylab
% - the title of this plot (titl)
% - the color (col) to be used to depict lines in this plot

% Invoked by: the user or MAIN.m
% Invokes: none

function [] = ave_cost_fig(gra,fig,betas,bbar,Jave,Jmin,xlab,ylab,titl,col)
loglog(betas,Jave','color',col,'linewidth',gra.lw)
hold on
grid on
yl = [10^floor(log10(Jmin)-1) 10^floor(log10(Jmin)+3)];
loglog([bbar bbar],[yl(1) Jmin],'k--','linewidth',gra.lw-0.5)
xl = xlim;
loglog([xl(1) bbar],[Jmin Jmin],'k--','linewidth',gra.lw-0.5)
plot([bbar bbar],[Jmin Jmin],'color',col,'marker','.',...
    'markeredgecolor',col*0.65,'markersize',50)
ylim(yl)


ax = fig;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = gra.ftsz;
title(titl,'fontsize',gra.ftsz,'interpreter','latex')
xlabel(xlab,'fontsize',gra.ftsz,'interpreter','latex')
ylabel(ylab,'fontsize',gra.ftsz,'interpreter','latex')

drawnow

end