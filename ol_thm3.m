% This function finds a solution for the proposed DDPC scheme in this 
% paper, which is based on the approximation discussed in Theorem 3; 
% it computes
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% - the associated control and output cost (prd.yu_cost)
% - the associated cost due to this regularization technique (prd.reg_cost)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given dataset contained in dpc
% - the prediction parameters and references contained in prd

% Invoked by: 
% - ol() in cl.m, to be started
% Invokes: none


function [prd] = ol_thm3(clx,dpc,prd,t)

Nprime = dpc.N-(clx.sys.m+clx.sys.p)*dpc.TTini;
prd.beta = dpc.N/Nprime; % adjustement: unbiased sample variance in (68)
prd.Reg_Thm3 = prd.beta*clx.omega;

if clx.opt.is_unconstrained
    opt = clx.opt;
    r1 = length(prd.gamma1);
    prd.gamma2 = dpc.L32'*opt.TQ*dpc.L32 + dpc.L22'*opt.TR*dpc.L22 + ...
        prd.Reg_Thm3(r1+1:end,r1+1:end);
    prd.gamma2 = (prd.gamma2^-1)*...
        (dpc.L32'*opt.TQ*(prd.yr-dpc.L31*prd.gamma1) + ...
        dpc.L22'*opt.TR*(prd.ur-dpc.L21*prd.gamma1)...
        -prd.Reg_Thm3(r1+1:end,1:r1)*prd.gamma1);
    prd = prd_uy(dpc,prd);
    
    prd.yu_cost = (prd.yr-prd.y)'*opt.TQ*(prd.yr-prd.y)+...
        (prd.ur-prd.u)'*opt.TR*(prd.ur-prd.u);
    prd.reg_cost = prd.beta*norm(prd.gamma2)^2;
else
    prd = cvx_sol(clx,dpc,prd,t,7);
end

end
