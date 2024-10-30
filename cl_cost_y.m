% This is an AUXILIARY function; it computes 
% - the closed-loop cost related to the output only
% Given
% - the optimization parameters and references contained in opt
% - the output y

% Invoked by: 
% - cl.m, after the end of the closed-loop
% Invokes: none

function cost = cl_cost_y(opt,y)

cost = norm(opt.yr(:,1:size(y,2))-y,"fro")^2;
    
end