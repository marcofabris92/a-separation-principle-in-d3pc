% This is a UTILITY function to draw output trajectories.
% It is used to look for "weird trajectories" to be shown in Figure 2.(a)
% and Figure 2.(b). The user will use the result stored in jj to set lines
% 74 and 75 of track_err_fig3.
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


function jj = weird_trajectories(gra,clx,NMC,s,xlab,ylab,titl,col)
% s = signal, r = reference   Tv x NMC

% dimension of s: d_s, Tv, NMC
d_s = size(s,1);
Tv = clx.Tv;
tt = (1:Tv);

r = clx.opt.yr(:,tt);
jj = [];
for j = 1:NMC
    found = 0;
    t = 0;
    while t < Tv && ~found
        t = t + 1;
        if abs(r(t)-s(1,t,j)) > 5
            found = 1;
        end
    end
    if found
        jj = [jj j];
    end
end


gra_pos = gra.pos;
for i = 1:d_s
    figure('position',gra_pos)
    plot(tt-1,r(i,:),'k','linewidth',gra.lw)
    hold on
    grid on
    
    
    for j = 1:length(jj)
        u = jj(j);
        plot(tt-1,s(1,:,u),'-','linewidth',1)
    end
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