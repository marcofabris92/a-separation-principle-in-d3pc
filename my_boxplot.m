% This is a UTILITY function to draw boxplots nicely.
% Given
% - graphic options gra
% - the distribution names x
% - the distribution values y
% - a vector of constants, just in case horizontal references are needed
%   to be depicted as well
% - a vector legs containing legend items (usually empty)
% - the legend location loc (usually empty)
% - the x label xlab (usually empty)
% - the y label ylab (usually empty)
% - the title of this plot (titl)
% - a flag Ylog (if Ylog = 1 then logarithmic scale on the y axis is
%   activated)
% - a flag mdn (if mdn = 1 then the median of all distribution vales y is
%   depicted as well)

% Invoked by: the user or MAIN.m
% Invokes: none


function fig = ...
    my_boxplot(gra,x,y,constants,legs,loc,xlab,ylab,titl,Ylog,mdn)

figure('position',gra.pos)
y = squeeze(y)';
boxplot(y,x)
hold on
grid on
if Ylog
    yl = ylim;
    ylim([0 yl(end)])
    ax = gca;
    ax.YAxis.Scale ='log';
end
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = gra.ftsz;
xl = xlim;
yl = ylim;
if length(x) > 10
    xtks = (1:5:x(end));
    set(gca, 'xtick',xtks,'XTickLabels',[0 5:5:x(end)])
end
gold = [255 215 0]/255;
h_const = [];
for j = 1:length(constants)
    sty = '-';
    color_map = {'g', 'm', gold, 'c', 'b'};
    col = [rand rand rand];
    if j > 5
        sty = '--';
    end
    if j > 10
        sty = '.-';
    end
    if j > 15
        sty = '*';
    else
        col = color_map{j};
    end
    c = constants(j);
    if c < yl(1)
        yl(1) = c;
    end
    if c > yl(end)
        yl(end) = c;
    end
    h = plot([xl(1) xl(end)],[c c],...
        'color',col,'linestyle',sty,'linewidth',gra.lw-0.5);
    h_const = [h_const h];
    if Ylog
        ax = gca;
        ax.YAxis.Scale ='log';
    end
end
if mdn
    plot(median(y),'mo-','linewidth',gra.lw-1)
end
ylim(yl)
xlabel(xlab,'fontsize',gra.ftsz,'interpreter','latex')
ylabel(ylab,'fontsize',gra.ftsz,'interpreter','latex')
title(titl,'fontsize',gra.ftsz,'interpreter','latex')
if ~isempty(legs)
    legend(h_const,legs,'fontsize',gra.ftsz,'interpreter','latex',...
        'location',loc)
end

fig = gcf;
drawnow

end