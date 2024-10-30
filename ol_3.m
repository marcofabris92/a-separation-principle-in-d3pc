% This function finds a solution for the gamma-DDPC scheme first 
% introduced in [9] (exploiting just the regularization on gamma3); 
% it computes: 
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given data set contained in dpc
% - the prediction parameters and references contained in prd
% - the hyperparameter beta = beta3, if beta3 has been fixed; beta = [],
%   if beta3 needs to be tuned with the strategy proposed in [9]
% - the discrete time instant t-1 (i.e. t-th iteration of the closed-loop)

% Invoked by: 
% - ol() in cl.m, to be started
% - tuning.m, to be run with a fixed beta3
% Invokes:
% - cvx_sol.m, if optimization constraints are enforced
% - prd_uy.m, to compute prd.u and prd.y
% - tuning.m, to start the tuning of beta3
% - z_star.m, to compute gamma2 and gamma3


function [prd] = ol_3(clx,dpc,prd,beta,t)

if isempty(beta)
    prd = tuning(clx,dpc,prd,3,t);
else
    prd.beta = beta;
end

if clx.opt.is_unconstrained
    opt = clx.opt;
    W_star = [clx.W_star clx.W12*[zeros(opt.mT,opt.pT); dpc.L33]];
    gamma23 = W_star'*W_star;
    if prd.beta == +Inf
        W_star = gamma23(1:opt.mT,1:opt.mT)^-1;
        gamma23 = zeros(opt.mT+opt.pT,1);
        gamma23(1:opt.mT) = W_star*clx.W_star'*z_star(clx,dpc,prd);
    elseif prd.beta == 0
        gamma23 = pinv(gamma23)*W_star'*z_star(clx,dpc,prd);
    else
        for s = (1:opt.pT)+opt.mT
            gamma23(s,s) = gamma23(s,s) + prd.beta;
        end
        gamma23 = (gamma23^-1)*W_star'*z_star(clx,dpc,prd);
    end
    prd.gamma2 = gamma23(1:opt.mT);
    prd.gamma3 = gamma23(opt.mT+1:end);    
    prd = prd_uy(dpc,prd);
else
    prd = cvx_sol(clx,dpc,prd,t,3);
end

end