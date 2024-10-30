% This function finds a solution for the gamma-DDPC scheme first 
% introduced in [7] (exploiting both the regularizing terms); 
% it computes: 
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given data set contained in dpc
% - the prediction parameters and references contained in prd
% - the hyperparameter beta = [beta2 beta3], if [beta2 beta3] have been 
%   fixed; beta = [], if [beta2 beta3] need to be tuned with the 
%   strategy proposed in [7]

% Invoked by: 
% - ol() in cl.m, to be started
% Invokes:
% - aug_diag(), see below
% - cvx_sol.m, if optimization constraints are enforced
% - prd_uy.m, to compute prd.u and prd.y
% - z_star.m, to compute gamma2 and gamma3


function [prd] = ol_23(clx,dpc,prd,beta,t)

if isempty(beta)
    error('Online tuning not coded yet for ol_23.')
else
    prd.beta = beta;
end

if clx.opt.is_unconstrained
    opt = clx.opt;
    beta2 = beta(1);
    beta3 = beta(2);
    
    % there are 9 cases below: (C3,C2), with C3 = 1,2,3 and C2 = 1,2,3
    % Ci = 1 --> beta_i is zero
    % Ci = 2 --> beta_i is finite and different from zero
    % Ci = 3 --> beta_i is infinite
    
    % (3,3) by default
    prd.gamma2 = zeros(opt.mT,1);
    prd.gamma3 = zeros(opt.pT,1);
    % (1,1)
    if beta2 == 0 && beta3 == 0
        prd.gamma2 = -(dpc.L22^-1)*dpc.L21*prd.gamma1;
        prd.gamma3 = clx.L33_inv*(prd.yr-...
            (dpc.L31*prd.gamma1+dpc.L32*prd.gamma2));
    end
    % (1,2)
    if beta2 > 0 && beta2 < +Inf && beta3 == 0
        Finv = aug_diag(dpc.F0,beta2)^-1;
        prd.gamma3 = pinv(dpc.H0-dpc.G'*Finv*dpc.G)*...
            (prd.c3-dpc.G'*Finv*prd.c2);
        prd.gamma2 = Finv*(prd.c2-dpc.G*prd.gamma3);
    end
    % (1,3)
    if beta2 == +Inf && beta3 == 0
        prd.gamma3 = pinv(dpc.H0)*prd.c3;
    end  
    % (2,1)
    if beta2 == 0 && beta3 > 0 && beta3 < +Inf
        H_inv = aug_diag(dpc.H0,beta3)^-1;
        prd.gamma2 = ((dpc.F0-dpc.G*H_inv*dpc.G')^-1)*...
            (prd.c2-dpc.G*H_inv*prd.c3);
        prd.gamma3 = H_inv*(prd.c3-dpc.G'*prd.gamma2);
    end
    % (2,2)
    if beta2 > 0 && beta2 < +Inf && beta3 > 0 && beta3 < +Inf
        H_inv = aug_diag(dpc.H0,beta3)^-1;
        prd.gamma2 = ((aug_diag(dpc.F0,beta2)-dpc.G*H_inv*dpc.G')^-1)*...
            (prd.c2-dpc.G*H_inv*prd.c3);
        prd.gamma3 = H_inv*(prd.c3-dpc.G'*prd.gamma2);
    end
    % (2,3)
    if beta2 == +Inf && beta3 > 0 && beta3 < +Inf
        prd.gamma3 = (aug_diag(dpc.H0,beta3)^-1)*prd.c3;
    end
    % (3,1)
    if beta2 == 0 && beta3 == +Inf
        prd.gamma2 = (dpc.F0^-1)*prd.c2;
    end
    % (3,2)
    if beta2 > 0 && beta2 < +Inf && beta3 == +Inf
        prd.gamma2 = (aug_diag(dpc.F0,beta2)^-1)*prd.c2;
    end
    prd = prd_uy(dpc,prd);
else 
    prd = cvx_sol(clx,dpc,prd,4);
end

end



% This function is used to add a constant increment beta to the diagonal 
% of a given matrix M
% Given
% - the matrix M
% - the increment beta

% Invoked by:
% - ol_23.m, whenever the diagonal of a matrix needs to be incremented
% Invokes:none


function M = aug_diag(M,beta)

for s = 1:size(M,1)
    M(s,s) = M(s,s) + beta;
end

end