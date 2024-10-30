% This subroutine runs inner Monte Carlo loops to explore performance 
% indexes of all the offline approaches;
% it was coded just to isolate such iterations from outerMC.m and tidy this
% whole code up
% Given 
% - several parameters, data and variables passed by outerMC.m

% Invoked by:
% - outerMC.m, to start a grid optimization in order to perform the offline 
%   tuning of certain methods, such as gamma-DDPC and DeePC
% Invokes:
% - cl.m, to start each closed-loop over which tuning parameters are tested
% - sys_noise.m, to return the j_inner-th noise realization


function [J1_grid,J2_grid,J3_grid,J23_grid,exe_time] = innerMC_grids...
    (sys,clx,Lpts,dpc,idx1,idx2,idx_1,idx_2,J1_grid,J2_grid,J3_grid,...
    J23_grid,exe_time,j)

% j_inner-th noise realization
[sys.vx,sys.vy] = sys_noise(sys,sys.TvT);
clx.sys = sys;


% i-th grid optimization for beta2 and beta3 (gamma-DDPC)
%fprintf('\ngrid optimization for beta2 and beta3\n')
tStart = tic;

parfor i = 1:Lpts % parfor
    if mod(i,1000) == 0 % just to show iterations (reduce that 1000)
        fprintf(['i = ' num2str(i) '\n'])
    end
    sol = cl(clx,dpc,2,0,clx.beta2_(i));
    J2_grid(i,j) = J2_grid(i,j) + sol.J;
    sol = cl(clx,dpc,3,0,clx.beta3_(i));
    J3_grid(i,j) = J3_grid(i,j) + sol.J;
end
exe_time.b2and3g = exe_time.b2and3g + toc(tStart);


% joint i1-i2-th grid optimization for beta2 and beta3 (gamma-DDPC)
%fprintf('\njoint grid optimization for beta23\n')

I2 = 1:length(idx2);
tStart = tic;
parfor i1 = 1:length(idx1) % parfor
    for i2 = I2
        sol = cl(clx,dpc,4,0,[clx.beta2_(i1) clx.beta3_(i2)]);
        J23_grid(i1,i2,j) = J23_grid(i1,i2,j) + sol.J;
    end
end
exe_time.b23g = exe_time.b23g + toc(tStart);


% joint i1-i2-th grid optimization for lambda1 and lambda2 (DeePC)
%fprintf('\ngrid optimization for DeePC\n')
I2 = 1:length(idx_2);
tStart = tic;
parfor i1 = 1:length(idx_1) % parfor
    for i2 = I2
        %[clx.lmb1_(i1) clx.lmb2_(i2)]
        sol = cl(clx,dpc,5,0,[clx.lmb1_(i1) clx.lmb2_(i2)]);
        J1_grid(i1,i2,j) = J1_grid(i1,i2,j) + sol.J;
    end
end
exe_time.l12g = exe_time.l12g + toc(tStart);
              
end