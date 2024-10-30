% This function is employed only within gamma-DDPC schemes [7], [8], [9] 
% in order to provide an intermediate result that is needed for the 
% subsequent computation of the variables gamma2 and gamma3
% Given
% - the closed-loop parameters and settings contained in clx
% - the parameters related to the given dataset contained in dpc
% - the prediction parameters and references contained in prd

% Invoked by:
% - ol_2
% - ol_3
% - ol_GT
% Invokes: none


function value = z_star(clx,dpc,prd)

value = clx.W12*([prd.ur; prd.yr]-[dpc.L21; dpc.L31]*prd.gamma1);

end