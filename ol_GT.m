% This is an AUXILIARY function that finds noiseless solutions for the
% gamma-DDPC schemes [7], [8], [9] in order to provide a ground truth (GT);
% it computes: 
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given data set contained in dpc
% - the prediction parameters and references contained in prd
% - the discrete time instant t-1 (i.e. t-th iteration of the closed-loop)

% Invoked by: 
% - ol() in cl.m, to be started (whenever the GT is optionally required)
% Invokes:
% - cvx_sol.m, if optimization constraints are enforced
% - prd_uy.m, to compute prd.u and prd.y
% - z_star.m, to compute gamma2


function [prd] = ol_GT(clx,dpc,prd,t)

if clx.opt.is_unconstrained
    prd.gamma2 = pinv(clx.W_star)*z_star(clx,dpc,prd);
    prd = prd_uy(dpc,prd); 
else
    prd = cvx_sol(clx,dpc,prd,t,0);
end

end