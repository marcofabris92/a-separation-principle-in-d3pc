% This function is an implementation of DeePC according to the notation
% used in [15]; it computes: 
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given data set contained in dpc
% - the prediction parameters and references contained in prd
% - the optimization parameters and references contained in opt

% Invoked by: 
% - ol() in cl.m, to be started
% Invokes:
% - cvx_sol.m, if optimization constraints are enforced
% - objective (see documentation below)


function [prd] = ol_DeePC(clx,dpc,prd,t)

if clx.opt.is_unconstrained

if prd.lmb1 == +Inf
    my_options = prd.options;
    my_options.Algorithm = 'sqp';
    x_star = fmincon(@(x)objective(x,clx.opt,dpc,prd),[prd.g0; prd.g0],...
            clx.Ac,clx.bc,[dpc.Zp  zeros(2*dpc.Tini,dpc.N);
                           dpc.IPi zeros(dpc.N)],...
                [prd.zini; zeros(dpc.N,1)],[],[],[],my_options);
else
    x_star = fmincon(@(x)objective(x,clx.opt,dpc,prd),[prd.g0; prd.g0],...
            clx.Ac,clx.bc,[dpc.Zp zeros(2*dpc.Tini,dpc.N)],prd.zini,...
            [],[],[],prd.options);
end
prd.u = dpc.Uf*x_star(1:dpc.N);
prd.y = dpc.Yf*x_star(1:dpc.N);

else
    % careful: not sure this works correctly (CVX code is still under 
    % development; certainly the case in which lambda1 = +Inf is not 
    % handled yet)
    prd = cvx_sol(clx,dpc,prd,5);
end

end



% This function defines the objective function fval (along with its 
% gradient gval; its Hessian is a constant and thus it is computed in cl.m,
% see variable HH) to be provided to fmincon inside ol_DeePC.m
% Given
% - the optimization meta-variable x
% - the optimization parameters and references contained in opt
% - the parameters related to the given data set contained in dpc
% - the prediction parameters and references contained in prd

% Invoked by:
% - ol_DeePC.m, to run the optimization having a norm1 in its cost
% Invokes: none


function [fval,gval] = objective(x,opt,dpc,prd)

g = x(1:dpc.N);

if prd.lmb1 == +Inf
    fval = ...
        (prd.ur-dpc.Uf*g)'*opt.TR*(prd.ur-dpc.Uf*g) + ...
        (prd.yr-dpc.Yf*g)'*opt.TQ*(prd.yr-dpc.Yf*g) + ...
        prd.lmb2*sum(x(1+dpc.N:end));
else
    fval = ...
        (prd.ur-dpc.Uf*g)'*opt.TR*(prd.ur-dpc.Uf*g) + ...
        (prd.yr-dpc.Yf*g)'*opt.TQ*(prd.yr-dpc.Yf*g) + ...
        prd.lmb1*g'*dpc.IPi2*g + prd.lmb2*sum(x(1+dpc.N:end));
end

if nargout > 1 % meaning: if gradient is required
    if prd.lmb1 == +Inf
        gval = [...
            dpc.Yf'*opt.TQ*(dpc.Yf*g-prd.yr) + ...
            dpc.Uf'*opt.TR*(dpc.Uf*g-prd.ur); 
            prd.lmb2*ones(dpc.N,1)];
    else
        gval = [...
            dpc.Yf'*opt.TQ*(dpc.Yf*g-prd.yr) + ...
            dpc.Uf'*opt.TR*(dpc.Uf*g-prd.ur) + ...
            prd.lmb1*dpc.IPi2*g; 
            prd.lmb2*ones(dpc.N,1)];
    end
end

end