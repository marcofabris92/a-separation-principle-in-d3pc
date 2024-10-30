% This is an AUXILIARY function; it computes 
% - the closed-loop cost related to the input only
% Given
% - the optimization parameters and references contained in opt
% - the input u

% Invoked by: 
% - cl.m, after the end of the closed-loop
% Invokes: none

function cost = cl_cost_u(opt,u)

cost = norm(opt.ur(:,1:size(u,2))-u,"fro")^2;
    
end