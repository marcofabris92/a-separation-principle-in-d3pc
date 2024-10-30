% This function finds a solution for the gamma-DDPC scheme first 
% introduced in [9]; it computes: 
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given data set contained in dpc
% - the prediction parameters and references contained in prd
% - the hyperparameter beta = beta2, if beta2 has been fixed; beta = [],
%   if beta2 needs to be tuned with the strategy proposed in [9]

% Invoked by: 
% - ol() in cl.m, to be started
% Invokes:
% - cvx_sol.m, if optimization constraints are enforced
% - prd_uy.m, to compute prd.u and prd.y
% - tuning.m, to start the tunung of beta2
% - z_star.m, to compute gamma2


function [prd] = ol_2(clx,dpc,prd,beta,t)

if isempty(beta)
    prd = tuning(clx,dpc,prd,2,t);
else
    prd.beta = beta;
end

if clx.opt.is_unconstrained
    opt = clx.opt;
    prd.gamma2 = clx.W_star'*clx.W_star;
    if prd.beta == +Inf
        prd.gamma2 = zeros(opt.mT,1);
    elseif prd.beta == 0
        prd.gamma2 = pinv(prd.gamma2)*clx.W_star'*z_star(clx,dpc,prd);
    else
        for s = 1:opt.mT
            prd.gamma2(s,s) = prd.gamma2(s,s) + prd.beta;
        end
        prd.gamma2 = (prd.gamma2^-1)*clx.W_star'*z_star(clx,dpc,prd);
    end
    prd = prd_uy(dpc,prd);
else 
    prd = cvx_sol(clx,dpc,prd,2);
end

end