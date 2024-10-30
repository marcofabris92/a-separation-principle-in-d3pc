% This is a UTILITY function to draw output trajectories.
% Version 2: used to get a more precise representation by selecting the 
% percentile value perc (see line 43) w.r.t. Version 1.
% Given
% - graphic options gra
% - the closed-loop variables, parameters and settings contained in clx
% - the total number of Monte Carlo runs
% - the distribution of a signal s to be depicted
% - the x label xlab
% - the y label ylab
% - the title of this plot (titl)
% - the color (col) to be used to depict lines in this plot

% Invoked by: the user or MAIN.m
% Invokes: none


function [] = track_err_fig2(gra,clx,NMC,s,xlab,ylab,titl,col)
% s = signal, r = reference   Tv x NMC

d_s = size(s,1);
Tv = clx.Tv;
if NMC == 1
    z = zeros(d_s,Tv,2);
    for i = 1:d_s
        for j = 1:Tv
            z(i,j,1) = s(i,j);
        end
    end
    s = z;
end

tt = (1:Tv);
r = clx.opt.yr(:,tt);
s_av = s(:,:,1);
if NMC > 1
    s_av = mean(s,Ls(s));
end

s_plus = zeros(d_s,Tv);
s_minus = zeros(d_s,Tv);
perc = 95;
j_star = round(perc/100*NMC);
if j_star == 0
    j_star = 1;
end
for t = tt
    tr_err = zeros(d_s,NMC);
    for j = 1:NMC
        for i = 1:d_s
            tr_err(i,j) = abs(r(i,t)-s(i,t,j));
        end
    end
    for i = 1:d_s
        tr_err(i,:) = sort(tr_err(i,:));
        s_plus(i,t) = r(i,t)+tr_err(i,j_star);
        s_minus(i,t) = r(i,t)-tr_err(i,j_star);
    end
end

gra_pos = gra.pos;
for i = 1:d_s
    figure('position',gra_pos)
    fill([tt-1, fliplr(tt-1)],[s_minus(i,:), fliplr(s_plus(i,:))],col,... 
        'FaceColor',col,'FaceAlpha',0.3,...
        'EdgeColor',col,'EdgeAlpha',0.3);
    hold on
    grid on
    plot(tt-1,s_minus(i,:),'color',col,'linewidth',gra.lw-1)
    plot(tt-1,s_plus(i,:),'color',col,'linewidth',gra.lw-1)
    plot(tt-1,r(i,:),'k','linewidth',gra.lw)
    %plot(tt-1,s(i,:,gra.j_star),'-','color',col/2,'linewidth',gra.lw)
    plot(tt-1,s_av(i,:),'--','color',col,'linewidth',gra.lw)
    xlim([1 Tv])
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = gra.ftsz;

    ylab_i = ylab;
    if d_s > 1
        ylab_i = [ylab ', component ' num2str(i)];
    end
    xlabel(xlab,'fontsize',gra.ftsz,'interpreter','latex')
    ylabel(ylab_i,'fontsize',gra.ftsz,'interpreter','latex')
    title(titl,'fontsize',gra.ftsz,'interpreter','latex')
end

drawnow

end
