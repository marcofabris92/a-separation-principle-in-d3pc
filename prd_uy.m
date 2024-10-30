% This function applies only to gamma-DDPC schemes [7], [8], [9]; 
% it computes:
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% Given
% - the parameters related to the given dataset contained in dpc
% - the prediction parameters and references contained in prd

% Invoked by: 
% - ol_2.m
% - ol_3.m
% - ol_23.m
% - ol_GT.m
% - ol_reg.m
% Invokes: none


function [prd] = prd_uy(dpc,prd)

prd.u = dpc.L21*prd.gamma1+dpc.L22*prd.gamma2;
prd.y = dpc.L31*prd.gamma1+dpc.L32*prd.gamma2;

end