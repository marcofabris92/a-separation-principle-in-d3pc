% This function finds a solution for the proposed DDPC scheme in this 
% paper, which is based on the Final Control Error; it computes
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% - the associated control and output cost (prd.yu_cost)
% - the associated cost due to this regularization technique (prd.reg_cost)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the prediction parameters and references contained in prd
% - the discrete time instant t-1 (i.e. t-th iteration of the closed-loop)

% Invoked by: 
% - ol() in cl.m, to be started
% Invokes: none


function [prd] = ol_FCE(clx,dpc,prd,t)

opt = clx.opt;
pbs = clx.pbs;

%% computing the next input to be applied
y_past = pbs.Xi_hat*prd.zini; 
y_constant = pbs.I_PsiY_hat*prd.yr-y_past;

if clx.opt.is_unconstrained
    prd.u = ((pbs.PsiU_hat'*pbs.Q*pbs.PsiU_hat+...
        opt.TR+clx.beta_0*pbs.Cuu)^-1)*...
        (pbs.PsiU_hat'*pbs.Q*y_constant+opt.TR*prd.ur...
        -clx.beta_0*(pbs.Cuz*prd.zini+pbs.Cuy*prd.yr));
    prd.y = pbs.I_PsiY_hat_inv*(y_past+pbs.PsiU_hat*prd.u);
else
    pbs.yc = y_constant;
    pbs.yp = y_past;
    clx.pbs = pbs;
    prd = cvx_sol(clx,dpc,prd,t,6);
end


%% returning values

prd.yu_cost = (y_constant-pbs.PsiU_hat*prd.u)'*pbs.Q*...
    (y_constant-pbs.PsiU_hat*prd.u)+(prd.ur-prd.u)'*opt.TR*(prd.ur-prd.u);
prd.reg_cost = clx.beta_0*[prd.zini; prd.u; prd.yr]'*pbs.LambdaDopt*...
    [prd.zini; prd.u; prd.yr];
    

end