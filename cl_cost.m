% This function computes 
% - the overall closed-loop cost
% Given
% - the optimization parameters and references contained in opt
% - the input u
% - the output y

% Invoked by: 
% - cl.m, after the end of the closed-loop
% Invokes: none

function cost = cl_cost(opt,u,y)

cost = 0;
Tv = size(u,2);
for t = 1:Tv
    ud_t = opt.ur(:,t)-u(:,t);
    yd_t = opt.yr(:,t)-y(:,t);
    cost = cost + ud_t'*opt.R*ud_t + yd_t'*opt.Q*yd_t;
end
    
end