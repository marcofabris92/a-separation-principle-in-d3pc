% This function computes the quantities reported in Proposition 3 of [9]
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the prediction parameters and references contained in prd

% Invoked by:
% - In tuning.m, to check whether |difference| is the smallest and use
%   num and den to find an appropriate value for best_beta
% Invokes: none


function [difference,num,den] = SD_3(clx,prd)

% num: ||\tilde{gamma3}||^2 was improved into
% num: ||(\hat{H}_s * kron(I_T,Lambda))^-1 * L33 * \tilde{gamma3}||^2
% (see data_generation.m and outerMC.m)
num = norm(clx.sys.HsLam_hat_inv_L33*prd.gamma3)^2;

% den: (T/N) * ||gamma12||^2
den = clx.TN*(norm([prd.gamma1; prd.gamma2])^2);

difference = num-den;

end