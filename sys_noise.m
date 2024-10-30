% This function is used to generate Gaussian noise vectors; it returns:
% - the output noise vy of length T
% - the process noise vx of length T
% Given
% - the statistics contained in sys
% - the lenght T

% Invoked by:
% - dynamics.m, to generate process and output noise for the model
% - innerMC_grids.m: to return the j_inner-th noise realization
% - outerMC.m, to return the noise realization of the j-th Monte Carlo run
% Invokes: none


function [vx,vy] = sys_noise(sys,T)

if sys.flag_noise
    vx = sys.V12*randn(size(sys.V12,2),T);
    vy = vx(sys.n+1:end,:);
    vx = vx(1:sys.n,:);
else
    vx = zeros(sys.n,T);
    vy = zeros(sys.p,T);
end

end