% This function performs tuning for the gamma-DDPC methods proposed in 
% [9]; it finds the best_beta(2 or 3) by minimizing |SD(2 or 3)| at each
% feedback step t
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given dataset contained in dpc
% - the prediction parameters and references contained in prd
% - the type of gamma-DDPC to be optimized (2 for beta2 or 3 for beta3)

% Invoked by:
% - ol_2.m, to perform tuning of beta2
% - ol_3.m, to perform tuning of beta3
% Invokes:
% - SD_2.m, to get the quantities in the approximation of (21) in [9]
% - SD_3.m, to get the quantities reported in Proposition 3 of [9]


function [prd] = tuning(clx,dpc,prd,type,t)

switch type
    case 2
        beta_ = clx.beta2_;
    case 3
        beta_ = clx.beta3_;
end

best_abs_diff = +Inf;
best_beta = NaN;
Lpts = length(beta_);
diffs = zeros(1,Lpts);
num = zeros(1,Lpts);
den = zeros(1,Lpts);
sign_num_R = 0;
sign_den_R = 0;
sign_num_L = 0;
sign_den_L = 0;
intersection_found = 0;
i_best = 0;
is_last_beta_finite = (beta_(end) < Inf);

% Main loop to best satisfy the propsed tuning strategies
% for beta2: see the approximation of (21) in [9]
% for beta3: see Proposition 3 of [9]
% The search for both betas is carried out in decreasing order
i = Lpts;
while ~intersection_found && i > 0
    % compute diffs, num and den
    switch type
        case 2
            [diffs(i),num(i),den(i)] =...
                SD_2(clx,ol_2(clx,dpc,prd,beta_(i),t));
        case 3
            [diffs(i),num(i),den(i)] =...
                SD_3(clx,ol_3(clx,dpc,prd,beta_(i),t));
    end
  
    % use num and den to find a suitable candidate for best_beta
    % adopted strategy: project the intersection of a cross formed by
    % num(i) den(i) num(i+1) den(i+1), which are quantities computed in
    % function of beta onto the beta axis
    if num(i) > den(i)
        sign_num_L = 1;
        sign_den_L = -1;
    elseif num(i) < den(i)
        sign_num_L = -1;
        sign_den_L = 1;
    else
        best_beta = beta_(i);
        intersection_found = 1;
    end
    
    if i < Lpts && ~intersection_found &&...
            sign_num_L*sign_num_R == -1 && sign_den_L*sign_den_R == -1
        del = [num(i+1)-den(i+1) num(i)-den(i)]';
        best_beta = ([beta_(i) -beta_(i+1)]*del)/([1 -1]*del); 
        intersection_found = 1;
    end
    
    if ~intersection_found && (i < Lpts || is_last_beta_finite)
        if num(i) > den(i)
            sign_num_R = 1;
            sign_den_R = -1;
        elseif num(i) < den(i)
            sign_num_R = -1;
            sign_den_R = 1;
        end
    end
   
    % check whether |diffs(i)| is the smallest over all i = 1:Lpts
    abs_diff_i = abs(diffs(i));
    if abs_diff_i < best_abs_diff
        best_abs_diff = abs_diff_i;
        i_best = i;
    end
    
    i = i - 1;
end

if ~intersection_found
    best_beta = beta_(i_best);
end

% filtering can be chosen to refine best_beta
if clx.LP 
    best_beta = ((t-1)*prd.beta+best_beta)/t;
end

prd.beta = best_beta;


% uncomment to see what is going on with the beta search at each step i
% figure
% semilogx(beta_,num,'k','linewidth',2)
% hold on
% grid on
% semilogx(beta_,den,'r','linewidth',1)
% title(['t = ' num2str(t) ', best beta' num2str(type) ' = ' ...
%     num2str(best_beta) ', i best = ' num2str(i_best)])
% pause
% close 

end