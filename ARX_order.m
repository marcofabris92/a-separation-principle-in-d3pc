% This function performs the ARX estimation and finds the best order
% according to the selected criterion; it returns
% - rho_star: the list of ARX orders
% - rho_min: the least ARX model order for this search (default: 1)
% - rho_max: the highest ARX model order for this search (default: Ndata-T)
% - nsv: an average of the estimated noise variances
% Given
% - the current version of the system parameters and references (sys)
% - the data sets stored in dat
% - the receding horizon T
% - the adopted criterion IC to validate ARX models
% - the labels, i.e. the set of indexes identifying each Monte Carlo run 
%   (it is possible to assign a vector of indexes to the variable labels)

% Invoked by:
% - rho_selection.m, to be launched
% Invokes: none


function [rho_star,rho_min,rho_max,nsv] = ARX_order(sys,dat,T,IC,labels)

% choose rho_min (>= 1)
rho_min = 1;

% choose rho_max (<= dat.Ndata-T)
rho_max = dat.Ndata-T;

% finds rho_star, i.e. the optimal lag to reconstruct the past trajectory
% rho_star selection is performed through the criterion specified in IC
LJs = length(labels);
rho_star = zeros(LJs,1);
noisevar_est = zeros(LJs,1);

% uncomment if you wish to see the process completion
% fprintf('Process completion:\n0%%\n')
% delta = floor(LJs/10);
warning off;
for j_ = 1:LJs
    % if mod(j_,delta) == 0
    %     fprintf([num2str(floor(100*j_/LJs)) '%%\n'])
    % end
    
    j = labels(j_);
    smallest_score = +Inf;
    rho = rho_min;
    N = 1;
    r = [];
    score = [];
    % try/check is used here to handle potential failures of ARX
    % estimations where too large model orders are imposed
    try
        while N > 0 && rho <= rho_max
            N = rho_max-rho+1;
            if j == 0
                YUdata = iddata(dat.y(:,1:N)',dat.u(:,1:N)',sys.Ts);
            else
                YUdata = iddata(dat.y_noisy(:,1:N,j)',...
                    dat.u_noisy(:,1:N,j)',sys.Ts);
            end
            arx_model = arx(YUdata,[rho*sys.na rho*sys.nb sys.nk]);
            if strcmp(IC,'FPE')
                score_j = fpe(arx_model);
            else
                score_j = aic(arx_model,IC);
            end
            if score_j < smallest_score
                noisevar_est = arx_model.noisevariance;
                smallest_score = score_j;
                rho_star(j_) = rho;
            end
            r = [r rho];
            score = [score score_j];
            rho = rho + 1;
        end
    catch
    end
    
% uncomment if you wish to see why rho_star(j_) was chosen
%     figure
%     grid on
%     hold on
%     plot(r,score,'k')
%     plot([rho_star(j_) rho_star(j_)],[min(score) max(score)],'r')
%     pause

end
warning on;

rho_star = round(mean(rho_star));
nsv = mean(noisevar_est);
% fprintf(['Estimated rho_star: ' num2str(rho_star) '\n'])

end