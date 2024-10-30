% This subroutine runs the j-th iteration of the Monte Carlo method;
% it was coded just to isolate such iterations from MAIN.m and tidy this
% whole code up
% Given 
% - several paramteres, data and variables passed by MAIN.m

% Invoked by:
% - MAIN.m, to run the j-th iteration of the Monte Carlo method
% Invokes:
% - cl.m, to start the closed-loop over which the control acts
% - dpc_ini.m, to reset dpc at any Monte Carlo experiment (if requested)
% - Hankel.m, to compute the Hankel matrices for the j-th data set
% - lqdec_part.m, to perform the LQ decomposition of the Hankel matrices
%   of the j-th Monte Carlo experiment for the gamma-DDPC approaches
% - HsLam_estimation.m to compute the quantity (Hs*Lambda^(1/2))^-1 that
%   is used for the online tuning of gamma-DDPC methods in [9] (see SD_2.m
%   and SD_3.m)
% - innerMC_grids.m, to start a grid optimization in order to perform the 
%   offline tuning of certain methods, such as gamma-DDPC and DeePC
% - pbsidopt_ini.m, to initialize the struct pbs needed for the
%   computations in ol_FCE.m
% - reg_weight.m, to compute the regularization weights used in ol_thm3.m
% - sys_noise.m, to generate the j-th noise realization (these realizations
%   are trajectories for both the state and the output; they are generated
%   in advance so that they can be stored to play around with them later,
%   when the execution of the whole script is over)


function [clx,sys,dpc,exe_time,...
    J1_grid,J2_grid,J3_grid,J23_grid,J1_ave,J2_ave,J3_ave,J23_ave,...
    beta2_bar,beta2_bar_,beta3_bar,beta23_bar,lmb12_bar,...
    J_MPC,J_MPC_u,J_MPC_y,u0,y0,y0_TR,...
    J2_bar,Ju2_bar,Jy2_bar,u2_bar,y2_bar,y2_bar_TR,y2_bar_KF,...
    n2_gamma12_2_bar,...
    J3_bar,Ju3_bar,Jy3_bar,u3_bar,y3_bar,y3_bar_TR,y3_bar_KF,...
    n2_gamma12_3_bar,gamma3_bar,...
    J23_bar,Ju23_bar,Jy23_bar,u23_bar,y23_bar,y23_bar_TR,y23_bar_KF,...
    n2_gamma12_23_bar,gamma3_23_bar,...
    J1_bar,Ju1_bar,Jy1_bar,u1_bar,y1_bar,y1_bar_TR,y1_bar_KF,...
    JFCE,JFCE_u,JFCE_y,u_FCE,y_FCE,y_FCE_TR,y_FCE_KF,beta_FCE,...
    Jthm3_tuned,Ju_thm3_tuned,Jy_thm3_tuned,u_thm3_tuned,...
    y_thm3_tuned,y_thm3_tuned_TR,...
    y_thm3_tuned_KF,n2_gamma12_thm3_tuned,beta_thm3_tuned,...
    J2_tuned,Ju2_tuned,Jy2_tuned,u2_tuned,y2_tuned,y2_tuned_TR,...
    y2_tuned_KF,n2_gamma12_2_tuned,beta2_tuned,SD2abs_actual,...
    J3_tuned,Ju3_tuned,Jy3_tuned,u3_tuned,y3_tuned,y3_tuned_TR,...
    y3_tuned_KF,n2_gamma12_3_tuned,gamma3_tuned,beta3_tuned,...
    SD3abs_actual...
    ] = outerMC(...
    clx,sys,dpc,...
    J1_grid,J2_grid,J3_grid,J23_grid,J1_ave,J2_ave,J3_ave,J23_ave,...
    beta2_bar,beta2_bar_,beta3_bar,beta23_bar,lmb12_bar,...
    J_MPC,J_MPC_u,J_MPC_y,u0,y0,y0_TR,...
    J2_bar,Ju2_bar,Jy2_bar,u2_bar,y2_bar,y2_bar_TR,y2_bar_KF,...
    n2_gamma12_2_bar,...
    J3_bar,Ju3_bar,Jy3_bar,u3_bar,y3_bar,y3_bar_TR,y3_bar_KF,...
    n2_gamma12_3_bar,gamma3_bar,...
    J23_bar,Ju23_bar,Jy23_bar,u23_bar,y23_bar,y23_bar_TR,y23_bar_KF,...
    n2_gamma12_23_bar,gamma3_23_bar,...
    J1_bar,Ju1_bar,Jy1_bar,u1_bar,y1_bar,y1_bar_TR,y1_bar_KF,...
    JFCE,JFCE_u,JFCE_y,u_FCE,y_FCE,y_FCE_TR,y_FCE_KF,beta_FCE,...
    Jthm3_tuned,Ju_thm3_tuned,Jy_thm3_tuned,u_thm3_tuned,...
    y_thm3_tuned,y_thm3_tuned_TR,...
    y_thm3_tuned_KF,n2_gamma12_thm3_tuned,beta_thm3_tuned,...
    J2_tuned,Ju2_tuned,Jy2_tuned,u2_tuned,y2_tuned,y2_tuned_TR,...
    y2_tuned_KF,n2_gamma12_2_tuned,beta2_tuned,SD2abs_actual,...
    J3_tuned,Ju3_tuned,Jy3_tuned,u3_tuned,y3_tuned,y3_tuned_TR,...
    y3_tuned_KF,n2_gamma12_3_tuned,gamma3_tuned,beta3_tuned,...
    SD3abs_actual,...
    dat,opt,gra,rho,idx1,idx2,idx_1,idx_2,rho_search_j,Lpts,NMC_inner,j)


visual = (mod(j,1) == 0);
if  gra.verbose && visual
    fprintf('\n__________________________\n')
    fprintf(['\nj = ' num2str(j) '\n'])
end

% resets dpc at each MC experiment if requested
if rho_search_j
    [dpc,sys] = dpc_ini(sys,dat,opt.T,rho(j));
    clx.Ac = [eye(dpc.N) -eye(dpc.N); -eye(dpc.N) -eye(dpc.N)];
    clx.bc = zeros(2*dpc.N,1);
    clx.TN = opt.T/dpc.N;
end

% j-th sigma2
clx.sigma2 = clx.nsvars(j);

% PBSIDopt preparation
tStart = tic;
clx.pbs = pbsidopt_ini(clx,dat,rho(j),j);
exe_time.trn_PBSIDopt = toc(tStart);

% noisy Hankel matrices
tStart = tic;
[dpc.Up,dpc.Uf] = Hankel(dpc,dat.u_noisy(:,:,j));
[dpc.Yp,dpc.Yf] = Hankel(dpc,dat.y_noisy(:,:,j));
dpc.H = [dpc.Up; dpc.Yp; dpc.Uf; dpc.Yf];
dpc.Zp = [dpc.Up; dpc.Yp];
dpc.Pi = pinv([dpc.Up; dpc.Yp; dpc.Uf])*[dpc.Up; dpc.Yp; dpc.Uf];
dpc.IPi = eye(dpc.N)-dpc.Pi;
dpc.IPi2 = dpc.IPi'*dpc.IPi;
exe_time.trn_DeePC = toc(tStart);

% LQ decomposition (noisy data)
dpc = lqdec_part(dpc,1);
dpc.F0 = dpc.L32'*opt.TQ*dpc.L32+dpc.L22'*opt.TR*dpc.L22;
dpc.G = dpc.L32'*opt.TQ*dpc.L33;
dpc.H0 = dpc.L33'*opt.TQ*dpc.L33;

% partial computation of repeatedly-used pieces of the solutions
clx.L33_inv = pinv(dpc.L33);
clx.W12 = blkdiag(opt.TR12,opt.TQ12);
clx.W_star = clx.W12*[dpc.L22; dpc.L32];
exe_time.trn_others = toc(tStart);

% estimation of (Hs*Lambda^(1/2))^-1 and related matrices
tStart = tic;
sys = HsLam_estimation(sys,dat,opt,rho(j),dat.Ndata,j);
sys.HsLam_hat_inv_L33 = sys.HsLam_hat_inv*dpc.L33;
exe_time.trn_23online = toc(tStart);
%check_Hs = norm(sys.Hs-sys.Hs_hat,"fro")

% preparation of the reg. weight (Thm3)
tStart = tic;
clx.omega = reg_weight(dpc,opt.T,opt.TQ);
exe_time.trn_S = toc(tStart);

    
%% grid optimization to perform offline tuning of regularization params.

exe_time.b2and3g = 0;
exe_time.b23g = 0;
exe_time.l12g = 0;
if ~sys.kind_of_setup
    if gra.verbose && visual
        fprintf('Optimization on grid ...\n')
    end
    for j_inner = 1:NMC_inner
        fprintf(['j_inner = ' num2str(j_inner) '\n'])
        [J1_grid,J2_grid,J3_grid,J23_grid,exe_time] = innerMC_grids(...
            sys,clx,Lpts,dpc,idx1,idx2,idx_1,idx_2,...
            J1_grid,J2_grid,J3_grid,J23_grid,exe_time,j);
    end

    exe_time.b2and3g = exe_time.b2and3g / NMC_inner;
    exe_time.b23g = exe_time.b23g / NMC_inner;
    exe_time.l12g = exe_time.l12g / NMC_inner;
    J1_grid(:,:,j) = J1_grid(:,:,j) / NMC_inner;
    J2_grid(:,j) = J2_grid(:,j) / NMC_inner;
    J3_grid(:,j) = J3_grid(:,j) / NMC_inner;
    J23_grid(:,:,j) = J23_grid(:,:,j) / NMC_inner;
    J1_ave = J1_ave + J1_grid(:,:,j);
    J2_ave = J2_ave + J2_grid(:,j);
    J3_ave = J3_ave + J3_grid(:,j);
    J23_ave = J23_ave + J23_grid(:,:,j);
    [il1_bar_j,il2_bar_j] = ... 
        (find(J1_grid(:,:,j) == min(min(J1_grid(:,:,j)))));
    [~,i2_bar_j] = min(J2_grid(:,j));
    [~,i3_bar_j] = min(J3_grid(:,j));
    [i22_bar_j,i33_bar_j] = ...
        (find(J23_grid(:,:,j) == min(min(J23_grid(:,:,j)))));
    beta2_bar(j) = clx.beta2_(i2_bar_j);
    beta2_bar_(j) = beta2_bar(j)/mean(diag(dpc.L33).^2);
    beta3_bar(j) = clx.beta3_(i3_bar_j);
    beta23_bar(:,j) = [clx.beta2_(i22_bar_j(1)) clx.beta3_(i33_bar_j(1))]';
    lmb12_bar(:,j) = [clx.lmb1_(il1_bar_j(1)) clx.lmb2_(il2_bar_j(1))]';
end

%% computing costs for the final boxplots
% j-th noise realization
[sys.vx,sys.vy] = sys_noise(sys,sys.TvT);
clx.sys = sys;

% KF-based oracle (MPC)
if  gra.verbose && visual
    fprintf('j-th KF-based oracle (MPC)\n')
end
tStart = tic;  
sol = cl(clx,dpc,1,0,[]);
exe_time.KF = toc(tStart);
J_MPC(j) = sol.J;
J_MPC_u(j) = sol.Ju;
J_MPC_y(j) = sol.Jy;
u0(:,:,j) = sol.uf;
y0(:,:,j) = sol.yf_hat;
y0_TR(:,:,j) = sol.yf;
if clx.gra.KF
    check = norm(sol.yf_hat-sol.yf,"fro");
    fprintf('************************************\n')
    fprintf('KALMAN-BASED NOISY ORACLE\n')
    fprintf(['Noisy oracle check: ' num2str(check) '\n'])
    fprintf(['sol.J = ' num2str(sol.J) '\n'])
    fprintf('************************************\n')
end
fprintf(['    JMPC_j = ' num2str(sol.J) '\n'])

% applying beta2 bar already found
fprintf('beta2 chosen a-posteriori\n')
tStart = tic;
sol = cl(clx,dpc,2,1,beta2_bar(j));
exe_time.b2 = toc(tStart);
J2_bar(j) = sol.J;
Ju2_bar(j) = sol.Ju;
Jy2_bar(j) = sol.Jy;
u2_bar(:,:,j) = sol.uf;
y2_bar(:,:,j) = sol.yf_hat;
y2_bar_TR(:,:,j) = sol.yf;
y2_bar_KF(:,:,j) = sol.yf_hat_star;
n2_gamma12_2_bar(:,j) = sol.n2gamma12;
fprintf(['    J2_j bar = ' num2str(sol.J) '\n'])

% applying beta3 bar already found
fprintf('beta3 chosen a-posteriori\n')
tStart = tic;
sol = cl(clx,dpc,3,1,beta3_bar(j));
exe_time.b3 = toc(tStart);
J3_bar(j) = sol.J;
Ju3_bar(j) = sol.Ju;
Jy3_bar(j) = sol.Jy;
u3_bar(:,:,j) = sol.uf;
y3_bar(:,:,j) = sol.yf_hat;
y3_bar_TR(:,:,j) = sol.yf;
y3_bar_KF(:,:,j) = sol.yf_hat_star;
n2_gamma12_3_bar(:,j) = sol.n2gamma12;
gamma3_bar(:,:,j) = sol.gamma3;
fprintf(['    J3_j bar = ' num2str(sol.J) '\n'])

% applying beta23 bar already found
fprintf('beta23 chosen a-posteriori\n')
tStart = tic;
sol = cl(clx,dpc,4,1,beta23_bar(:,j)');
exe_time.b23 = toc(tStart);
J23_bar(j) = sol.J;
Ju23_bar(j) = sol.Ju;
Jy23_bar(j) = sol.Jy;
u23_bar(:,:,j) = sol.uf;
y23_bar(:,:,j) = sol.yf_hat;
y23_bar_TR(:,:,j) = sol.yf;
y23_bar_KF(:,:,j) = sol.yf_hat_star;
n2_gamma12_23_bar(:,j) = sol.n2gamma12;
gamma3_23_bar(:,:,j) = sol.gamma3;
fprintf(['    J23_j bar = ' num2str(sol.J) '\n'])

% applying lambda12 bar already found
fprintf('lambda12 chosen a-posteriori\n')
tStart = tic;
sol = cl(clx,dpc,5,1,lmb12_bar(:,j)');
exe_time.l12 = toc(tStart);
J1_bar(j) = sol.J;
Ju1_bar(j) = sol.Ju;
Jy1_bar(j) = sol.Jy;
u1_bar(:,:,j) = sol.uf;
y1_bar(:,:,j) = sol.yf_hat;
y1_bar_TR(:,:,j) = sol.yf;
y1_bar_KF(:,:,j) = sol.yf_hat_star;
fprintf(['    J1_j bar = ' num2str(sol.J) '\n'])

% FCE tuning strategy
fprintf('Tuning through PBSIDopt\n')
tStart = tic;
sol = cl(clx,dpc,6,1,[]);
exe_time.O = toc(tStart);
JFCE(j) = sol.J;
JFCE_u(j) = sol.Ju;
JFCE_y(j) = sol.Jy;
u_FCE(:,:,j) = sol.uf;
y_FCE(:,:,j) = sol.yf_hat;
y_FCE_TR(:,:,j) = sol.yf;
y_FCE_KF(:,:,j) = sol.yf_hat_star;
beta_FCE(:,j) = sol.beta_FCE;
fprintf(['    JFCE_j tuned = ' num2str(sol.J) '\n'])

% best regularization according to Theorem 3
fprintf('Tuning for best regularization\n')
tStart = tic;
sol = cl(clx,dpc,7,1,[]);
exe_time.S = toc(tStart);
Jthm3_tuned(j) = sol.J;
Ju_thm3_tuned(j) = sol.Ju;
Jy_thm3_tuned(j) = sol.Jy;
u_thm3_tuned(:,:,j) = sol.uf;
y_thm3_tuned(:,:,j) = sol.yf_hat;
y_thm3_tuned_TR(:,:,j) = sol.yf;
y_thm3_tuned_KF(:,:,j) = sol.yf_hat_star;
n2_gamma12_thm3_tuned(:,j) = sol.n2gamma12;
beta_thm3_tuned(:,j) = sol.beta_thm3;
fprintf(['    Jthm3_j tuned = ' num2str(sol.J) '\n'])

% beta2 tuning strategy
fprintf('Tuning of beta2\n')
tStart = tic;
sol = cl(clx,dpc,2,1,[]);
exe_time.b2t = toc(tStart);
J2_tuned(j) = sol.J;
Ju2_tuned(j) = sol.Ju;
Jy2_tuned(j) = sol.Jy;
u2_tuned(:,:,j) = sol.uf;
y2_tuned(:,:,j) = sol.yf_hat;
y2_tuned_TR(:,:,j) = sol.yf;
y2_tuned_KF(:,:,j) = sol.yf_hat_star;
n2_gamma12_2_tuned(:,j) = sol.n2gamma12;
beta2_tuned(:,j) = sol.beta2_tuned;
SD2abs_actual(:,j) = abs(sol.actual_SD2);
fprintf(['    J2_j tuned = ' num2str(sol.J) '\n'])

% beta3 tuning strategy
fprintf('Tuning of beta3\n')
tStart = tic;
sol = cl(clx,dpc,3,1,[]);
exe_time.b3t = toc(tStart);
J3_tuned(j) = sol.J;
Ju3_tuned(j) = sol.Ju;
Jy3_tuned(j) = sol.Jy;
u3_tuned(:,:,j) = sol.uf;
y3_tuned(:,:,j) = sol.yf_hat;
y3_tuned_TR(:,:,j) = sol.yf;
y3_tuned_KF(:,:,j) = sol.yf_hat_star;
n2_gamma12_3_tuned(:,j) = sol.n2gamma12;
gamma3_tuned(:,:,j) = sol.gamma3;
beta3_tuned(:,j) = sol.beta3_tuned;
SD3abs_actual(:,j) = abs(sol.actual_SD3);
fprintf(['    J3_j tuned = ' num2str(sol.J) '\n'])

end