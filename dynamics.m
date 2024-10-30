% This function is used to compute updates of the underlying LTI system;
% it returns:
% - a windowed history of the output y over the current instant t and t+T-1
% - a windowed history of the state x over the current instant t and t+T
% - a copy of the generated/applied output noise vy 
% - a copy of the generated/applied process noise vx
% Given
% - the underlying LTI system parameters contained in sys
% - an initial condition xini at the current time t
% - the input u to be applied from t to t+T-1
% - the process noise vx to be applied from t to t+T (if empty, vx is
%   computed by sys_noise when it is required)
% - the output noise vy to be applied from t to t+T-1 (if empty, vy is
%   computed by sys_noise when it is required)

% Invoked by:
% - cl.m, while executing the closed-loop to be controlled
% - data_generation.m, to generate synthetic data and adjust parameters
% - dpc_ini.m, to (re)initalize the dpc struct when requested
% Invokes:
% - sys_noise.m, to generate the noise vectors vx and vy when requested

function [y,x,vy,vx] = dynamics(sys,xini,u,vx,vy)

T = size(u,2);
x = zeros(sys.n,T+1);
x(:,1) = xini;
y = zeros(sys.p,T);
if isempty(vx) && ~isempty(vy)
    vx = sys_noise(sys,T);
end
if ~isempty(vx) && isempty(vy)
    [~,vy] = sys_noise(sys,T);
end
if isempty(vx) && isempty(vy)
    [vx,vy] = sys_noise(sys,T);
end
for t = 1:T
    y(:,t) = sys.C*x(:,t)+sys.D*u(:,t)+sys.flag_noise*vy(:,t);
    x(:,t+1) = sys.A*x(:,t)+sys.B*u(:,t)+sys.flag_noise*vx(:,t);
end

end