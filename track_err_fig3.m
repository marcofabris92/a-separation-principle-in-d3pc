% This is a UTILITY function to draw output trajectories.
% Version 3: used to depict more specifically Figure 2 in this paper.
% To do so, you have to look for "weird trajectories" first by launching
% the function weird_trajectories.m and select manually jj at rows 
% 74 and 75 of this piece of code.
% Then, in order to draw Figures 2.(a) and 2.(b) activate lines 131 - 133 
% and deactivate line 134;
% instead, to draw Figures 2.(c) and 2.(d) activate line 134, 
% deactivate lines 131 - 133 and assign tt as you prefer at line 39.
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


function [] = track_err_fig3(gra,clx,NMC,s,xlab,ylab,titl,col)
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
q1 = round(25/100*NMC);
q3 = round(75/100*NMC);
w = 1.5;
if q1 == 0
    q1 = 1;
end
if q3 == 0
    q3 = 1;
end

jj = [];
for j = 1:NMC
    found = 0;
    t = 0;
    while t < Tv && ~found
        t = t + 1;
        if abs(r(t)-s(1,t,j)) > 1.401 % DeePC min_err = 5; 
            found = 1;                % FCE min_err = 1.401
        end
    end
    if found
        jj = [jj j];
    end
end
% weird trajectories - indexes for DeePC: 15    39    57
% weird trajectories - indexes for   FCE:  2    28    96
jj = [15 39 57];
jj = [2 28 96];

for t = tt
    tr_err = zeros(d_s,NMC);
    for j = 1:NMC
        for i = 1:d_s
            tr_err(i,j) = r(i,t)-s(i,t,j);
        end
    end
    for i = 1:d_s
        tr_err(i,:) = sort(tr_err(i,:));
        val_q1 = tr_err(i,q1);
        val_q3 = tr_err(i,q3);
        Wis_sup_can_val = val_q3+w*(val_q3-val_q1);
        Wis_inf_can_val = val_q1-w*(val_q3-val_q1);
        j = NMC;
        Wis_sup = 0;
        while Wis_sup == 0 && j <= NMC && j >= q3
            if tr_err(i,j) < Wis_sup_can_val
                Wis_sup = j;
            else
                j = j - 1;
            end
        end
        if j == q3 && Wis_sup == 0
            Wis_sup = q3;
        end
        j = 1;
        Wis_inf = 0;
        while Wis_inf == 0 && j >= 1 && j <= q1
            if tr_err(i,j) > Wis_inf_can_val
                Wis_inf = j;
            else
                j = j + 1;
            end
        end
        if j == q1 && Wis_inf == 0
            Wis_inf = q1;
        end
        s_plus(i,t) = r(i,t)-tr_err(i,Wis_sup);
        s_minus(i,t) = r(i,t)-tr_err(i,Wis_inf);
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
    plot(tt-1,r(i,:),'k','linewidth',gra.lw)
    %plot(tt-1,s(i,:,jj(1)),'-','color',[col/2 0.5],'linewidth',1)
    %plot(tt-1,s(i,:,jj(2)),'-','color',[col/2 0.5],'linewidth',1)
    %plot(tt-1,s(i,:,jj(3)),'-','color',[col/2 0.5],'linewidth',1)
    plot(tt-1,s(i,:,gra.j_star),'--','color',col/2,'linewidth',gra.lw)
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