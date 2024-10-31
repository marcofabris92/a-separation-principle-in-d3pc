% This script is used to execute the Numerical Examples conducted in 
% Section 7 of the paper entitled:

% "Harnessing Uncertainty for a Separation Principle in Direct
% Data-Driven Predictive Control"

% Authors:
% Prof. Alessandro Chiuso* Univ. of Padua      <alessandro.chiuso@unipd.it>
% Dr.   Marco Fabris**     Univ. of Padua         <marco.fabris.1@unipd.it>
% Prof. Valentina Breschi  Eindhoven Univ. of Technology <v.breschi@tue.nl>
% Prof. Simone Formentin   Polit. di Milano    <simone.formentin@polimi.it>
%
% *  A. Chiuso is the main algorithm developer of these Matlab scripts.
% ** M. Fabris is the main software developer of these MatLab scripts.

% Research Article submitted to the journal "Automatica", ID: 23-1641
% Publication history:
% - first submission: December 23rd, 2023
% - second submission: August 2nd, 2024
% - accepted on: October 3rd, 2024
% - last submission: October 31st, 2024

% Abstract:
% *************************************************************************
% Model Predictive Control (MPC) is a powerful method for complex system 
% regulation, but its reliance on an accurate model poses many limitations 
% in real-world applications. Data-driven predictive control (DDPC) aims at
% overcoming this limitation, by relying on historical data to provide 
% information on the plant to be controlled. 
% In this work, we present a unified stochastic framework for direct DDPC, 
% where control actions are obtained by optimizing the Final Control Error 
% (FCE), which is directly computed from available data only and 
% automatically weighs the impact of uncertainty on the control objective.
% Our framework allows us to establish a separation principle for 
% Predictive Control, elucidating the role that predictive models and their
% uncertainty play in DDPC. Moreover, it generalizes existing DDPC methods,
% like regularized Data-enabled Predictive Control (DeePC) and 
% $\gamma$-DDPC, providing a path toward noise-tolerant data-based control 
% with rigorous optimality guarantees. The theoretical investigation is 
% complemented by a series of numerical case studies, revealing that the 
% proposed method consistently outperforms or, at worst, matches existing 
% techniques without requiring tuning regularization parameters as other 
% methods do. 
% *************************************************************************

% This script and the remaining functions are documented at the beginning
% of each code fragment. See also the file README.txt for more details.

% Invoked by: the user
% Invokes:
% - ave_cost_fig.m (makes plots, see its documentation)
% - data_generatio.m, to generate data sets
% - dpc_ini.m, to initialize the dpc struct for the noise-free solution
% - get_GT.m, computes a noise-free solution to provide a ground truth (GT)
% - getPercIndex(), to find the index of a trajectory obtained from
%   a Monte Carlo experiment, given its cost and a percentile to select it
% - Hankel.m, to construct Hankel matrices for the noise-free solution
% - lqdec_part.m, LQ decomposition for the noise-free solution
% - model_selction.m, to select the model to be test
% - my_boxplot.m (makes plots, see its documentation)
% - outerMC.m, which is the body of the j-th iteration of the main Monte
%   Carlo experiment being conducted here
% - rho_selection.m, to select the model order for each Monte Carlo run
% - track_err_fig(1,2,3).m (makes plots, see its documentation)
% - weird_trajectories.m (makes plots, see its documentation)


clear all
close all
clc

%% execution time measurements
sys.kind_of_setup = 0; % 0 = setup 1 and 2; 1 = setup 3
global main_time cvx_sol_time cvx_sol_counter
cvx_sol_time = 0;
cvx_sol_counter = 0;
main_time = tic;

%% graphic options
savefig_flag = 1;
gra.lw = 2;
gra.ftsz = 19;
gra.verbose = 1;
gra.showIteratedFigs = 1;
gra.KF = 0;
nfigs = 1;
gra.show_ys_beta3_tuned = 0;
fig_counter = 1-nfigs;

%% activation of CVX for more complicated cases 
%  CVX was not used in our simulations, but it can be useful in general
% run ../cvx_startup.m (where .. has to be replaced with the correct path)
% cvx_solver sedumi % SDPT3
% cvx_precision best % low % medium % default % high % best

%% setup

% Structs:
% dat - contains info about data being used
% sys - contains info about the system model
% opt - contains info about the optimization parameters and cost functions
% clx - contains generic variables that need to be transferred to other
%    MATLAB functions during the optimization process

NMC = 1e2;                      % # of Monte Carlo samples (outer)
NMC_inner = 1e1;                % # of Monte Carlo samples (inner)
dat.Ndata = 250;             	% input/output series length (original 250)
sys.flag_noise = 1;           	% 1: sets noise in the model
        
opt.Wu = 5e-6;                  % input weight for cost (original 1e-2)
opt.Wy = 1;                     % output weight for cost (original 2*1e3)
% sys.q                         % state variance (noise in the model)
sys.r = 1*opt.Wu;               % output variance (noise in the model)
opt.T = 20;                    	% prediction horizon (original 20)
opt.I_T = eye(opt.T);           % identity matrix of dimension T
opt.O_T = zeros(opt.T);         % zero matrix of dimension T

opt.set_cnstr = '';             % list of set constraints for uf and yf_hat
opt.is_unconstrained = isempty(opt.set_cnstr);

% beta
beta_0s = [0.01 0.025 0.05 0.1 0.25 0.5 1 2.5 5 10];
clx.beta_0 = 1;                 % adjustment for FCE regularization
if clx.beta_0 == +Inf
    error("The case clx.beta_0 == +Inf has not been solved yet\n")
end

% generation of beta2 and beta3 grid
Lpts = 2*10^2;                      % # of points in the beta grids
clx.beta2_ = logspace(-3,1,Lpts);   % points for beta2
clx.beta3_ = logspace(-7,-3,Lpts);  % points for beta3 
Lpts = Lpts + 2;
clx.beta2_ = [0 clx.beta2_ +Inf];  
clx.beta3_ = [0 clx.beta3_ +Inf];


% decimation of indexes for betas 
% (unlock code below if you wish to reduce the number of points)
Lpts_r = Lpts;
idx1 = 1:Lpts;
idx2 = 1:Lpts;
% Lpts_r = min(clx.beta2_,clx.beta3_);
% ratio_Lpts = floor(Lpts/Lpts_r);
% i1 = 1; % <- select the step for beta2
% i2 = 1; % <- select the step for beta3
% idx1 = 1;
% idx2 = 1;
% for i1_ = 1:Lpts_r-1
%     i1 = i1 + ratio_Lpts;
%     i2 = i2 + ratio_Lpts;
%     idx1 = [idx1 i1];
%     idx2 = [idx2 i1];
% end


% generation of lambda1 and lambda2 grid
Lpts_r_lmb = 10;                         % # of points in the lambda grid
clx.lmb1_ = [0 logspace(-6,-1,Lpts_r_lmb-1) +Inf]; % points for lambda1
clx.lmb2_ = [0 logspace(-7,-2.5,Lpts_r_lmb)];      % points for lambda2
Lpts_r_lmb = Lpts_r_lmb + 1;                 


% decimation of indexes for lambdas
% (unlock code below if you wish to reduce the number of points)
idx_1 = 1:Lpts_r_lmb;
idx_2 = 1:Lpts_r_lmb;
% ratio_Lpts = floor(Lpts_lmb/Lpts_r_lmb);
% i1 = 1; % <- select the step for lambda1
% i2 = 1; % <- select the step for lambda2
% idx_1 = 1;
% idx_2 = 1;
% for i1_ = 1:Lpts_r_lmb-1
%     i1 = i1 + ratio_Lpts;
%     i2 = i2 + ratio_Lpts;
%     idx_1 = [idx_1 i1];
%     idx_2 = [idx_2 i1];
% end


clx.LP = 0;                     % low pass tuning heuristics
clx.gra = gra;                  % saving graphics options                  

sigma_u = 1;                    % std input training data
rho_search = 0;                 % 1: sets rho_hat via AICc once for all
rho_search_j = 1;               % 1: sets rho_hat via AICc at each
                                %    Monte Carlo experiment
                                
                                
% models (see model_selection.m)
models = [1]; % 0 and 1 have been tested in the proposed simulations
clr = 0; % to allow for filtering during data generation 
         % (see data_generation.m)

% information criteria: 'AICc', 'nAIC', 'AIC', 'BIC', 'FPE'
ICs = {'AICc'}; % <- choose your favorite criterias to determine rho_hat
exe_times = cell(NMC);

%% load setup and data
if sys.kind_of_setup
    location = ''; % <-- add your path
    filename = 'S0_model1_SNR20_wrkspc.mat'; % <-- choose the correct name
    if isempty(location)
        load(filename)
    else
        load([location filename]);
    end
    sys.kind_of_setup = 1;
end

clx.Tv = 500;    % Tv > 1              	


%% model loop
for model = models
 
    %% model selection:
    [sys,opt] = model_selection(sys,clx,opt,model); 
    clx.opt = opt;
   
    % SNR Var(y)/Var(noise) (in dB): list to be tried
    switch model
        case 0
            SNRs_dB = [20];
        case 1
            SNRs_dB = [20];
            clr = 1;
    end
    
    
for snr_target_dB = SNRs_dB
    fig_counter = fig_counter + nfigs;
        
    fprintf('////////////////////////////////////\n')
    fprintf(['Selected SNR: ' num2str(snr_target_dB) ' dB\n\n'])

    %% data generation
    % 1) adjusting sys.Q, sys.R, sys.S; computing Kalman paramters 
    % 2) initializing dat.u and dat.y (noiseless dataset)
    % 3) getting training dataset (NMC noisy datasets)
    [sys,dat] = data_generation(sys,dat,opt,sigma_u,snr_target_dB,clr,NMC);
    clx.opt = opt;
    
        
for ic = 1:length(ICs)
    
    % selection of the length of the past window (rho)
    IC = ICs{ic}; % chosen information criterion
    fprintf(['Estimating rho (' IC ') ...\n'])
    tStart = tic;
    [rho_GT,rho,sys,nsvars] = ...
        rho_selection(sys,dat,opt.T,IC,rho_search_j,rho_search,NMC);
    trn_rho = toc(tStart);
    fprintf('Values of rho found:')
    rho_GT
    rho         % contains the used trunctions for the model orders
    nsvars      % estimations of the noise variance
    average_rho = mean(rho);
    average_rho
    
    
    % Noise-free ground-truth solution
    % (not shown in the simulations; it can be used to understand limit
    % cases representing noise-free scenarios)
    copy_flag_noise = sys.flag_noise;
    sys.flag_noise = 0;
    [dpc,sys] = dpc_ini(sys,dat,opt.T,rho_GT);
    sys.flag_noise = copy_flag_noise;
    [dpc.UpGT,dpc.UfGT] = Hankel(dpc,dat.u);
    [dpc.YpGT,dpc.YfGT] = Hankel(dpc,dat.y); 
    dpc = lqdec_part(dpc,0);  
    clx.TN = opt.T/dpc.N;
    clx.sys = sys;
    sol_GT = get_GT(clx,dpc);
    
    
    % outer MC loop initialization
    % subfix 'bar' means 'minimum of the average'
    % subfix 'ave' means 'average of the minima'
    % J's: performance index
    % u's: stored input sequences
    % y's: stored output sequences
    % TR subfix: 
    % KF subfix:
    
    % oracle (a = MPC)
    J_MPC = zeros(1,NMC);
    J_MPC_u = zeros(1,NMC);
    J_MPC_y = zeros(1,NMC);
    
    u0 = zeros(opt.mT,clx.Tv,NMC);
    y0 = zeros(opt.pT,clx.Tv,NMC);
    y0_TR = zeros(opt.pT,clx.Tv,NMC);
    
    % DeePC method (a = DeePC)
    J1_grid = zeros(Lpts_r_lmb,Lpts_r_lmb,NMC);
    J1_ave = zeros(Lpts_r_lmb,Lpts_r_lmb);
    J1_bar = zeros(1,NMC);
    Ju1_bar = zeros(1,NMC);
    Jy1_bar = zeros(1,NMC);
    lmb12_bar = zeros(2,NMC);
    
    u1_bar = zeros(opt.mT,clx.Tv,NMC);
    y1_bar = zeros(opt.pT,clx.Tv,NMC);
    y1_bar_TR = zeros(opt.pT,clx.Tv,NMC);
    y1_bar_KF = zeros(opt.pT,clx.Tv,NMC);
    
    % gamma-DDPC method (a = gamma2, bar) + grid optimization
    J2_grid = zeros(Lpts,NMC);
    J2_ave = zeros(Lpts,1);
    J2_bar = zeros(1,NMC);
    Ju2_bar = zeros(1,NMC);
    Jy2_bar = zeros(1,NMC);
    beta2_bar = zeros(1,NMC);
    beta2_bar_ = zeros(1,NMC);
    
    u2_bar = zeros(opt.mT,clx.Tv,NMC);
    y2_bar = zeros(opt.pT,clx.Tv,NMC);
    y2_bar_TR = zeros(opt.pT,clx.Tv,NMC);
    y2_bar_KF = zeros(opt.pT,clx.Tv,NMC);
    n2_gamma12_2_bar = zeros(clx.Tv,NMC);
    
    % gamma-DDPC method (a = gamma3, bar) + grid optimization
    J3_grid = zeros(Lpts,NMC);
    J3_ave = zeros(Lpts,1);
    J3_bar = zeros(1,NMC);
    Ju3_bar = zeros(1,NMC);
    Jy3_bar = zeros(1,NMC);
    beta3_bar = zeros(1,NMC);
    
    u3_bar = zeros(opt.mT,clx.Tv,NMC);
    y3_bar = zeros(opt.pT,clx.Tv,NMC);
    y3_bar_TR = zeros(opt.pT,clx.Tv,NMC);
    y3_bar_KF = zeros(opt.pT,clx.Tv,NMC);
    n2_gamma12_3_bar = zeros(clx.Tv,NMC);
    gamma3_bar = zeros(opt.pT,clx.Tv,NMC);
    
    % gamma-DDPC method (a = gamma23, bar) + grid optimization
    J23_grid = zeros(Lpts_r,Lpts_r,NMC);
    J23_ave = zeros(Lpts_r,Lpts_r);
    J23_bar = zeros(1,NMC);
    Ju23_bar = zeros(1,NMC);
    Jy23_bar = zeros(1,NMC);
    beta23_bar = zeros(2,NMC);
    
    u23_bar = zeros(opt.mT,clx.Tv,NMC);
    y23_bar = zeros(opt.pT,clx.Tv,NMC);
    y23_bar_TR = zeros(opt.pT,clx.Tv,NMC);
    y23_bar_KF = zeros(opt.pT,clx.Tv,NMC);
    n2_gamma12_23_bar = zeros(clx.Tv,NMC);
    gamma3_23_bar = zeros(opt.pT,clx.Tv,NMC);
    
    % FCE method (a = FCE)
    JFCE_u = zeros(1,NMC);
    JFCE_y = zeros(1,NMC);
    JFCE = zeros(1,NMC);   
    
    beta_FCE = zeros(clx.Tv,NMC);
    
    u_FCE = zeros(opt.mT,clx.Tv,NMC);
    y_FCE = zeros(opt.pT,clx.Tv,NMC);
    yFCE_TR = zeros(opt.pT,clx.Tv,NMC);
    yFCE_KF = zeros(opt.pT,clx.Tv,NMC);
    n2_gamma12_o = zeros(clx.Tv,NMC);
    
    % approximated FCE regularizer (a = thm3)
    Ju_thm3_tuned = zeros(1,NMC);
    Jy_thm3_tuned = zeros(1,NMC);
    Jthm3_tuned = zeros(1,NMC);  
    
    beta_thm3_tuned = zeros(clx.Tv,NMC);
    
    u_thm3_tuned = zeros(opt.mT,clx.Tv,NMC);
    y_thm3_tuned = zeros(opt.pT,clx.Tv,NMC);
    y_thm3_tuned_TR = zeros(opt.pT,clx.Tv,NMC);
    y_thm3_tuned_KF = zeros(opt.pT,clx.Tv,NMC);
    n2_gamma12_thm3_tuned = zeros(clx.Tv,NMC);
    
    % gamma-DDPC method (a = gamma2, hat) + tuning strategy
    Ju2_tuned = zeros(1,NMC);
    Jy2_tuned = zeros(1,NMC);
    J2_tuned = zeros(1,NMC);   
    
    beta2_tuned = zeros(clx.Tv,NMC);
    SD2abs_actual = zeros(clx.Tv,NMC);
    
    u2_tuned = zeros(opt.mT,clx.Tv,NMC);
    y2_tuned = zeros(opt.pT,clx.Tv,NMC);
    y2_tuned_TR = zeros(opt.pT,clx.Tv,NMC);
    y2_tuned_KF = zeros(opt.pT,clx.Tv,NMC);
    n2_gamma12_2_tuned = zeros(clx.Tv,NMC);
    
    % gamma-DDPC method (a = gamma3, hat) + tuning strategy
    Ju3_tuned = zeros(1,NMC);
    Jy3_tuned = zeros(1,NMC);
    J3_tuned = zeros(1,NMC);    
    
    beta3_tuned = zeros(clx.Tv,NMC);
    SD3abs_actual = zeros(clx.Tv,NMC);
    
    u3_tuned = zeros(opt.mT,clx.Tv,NMC);
    y3_tuned = zeros(opt.pT,clx.Tv,NMC);
    y3_tuned_TR = zeros(opt.pT,clx.Tv,NMC);
    y3_tuned_KF = zeros(opt.pT,clx.Tv,NMC);
    n2_gamma12_3_tuned = zeros(clx.Tv,NMC);
    gamma3_tuned = zeros(opt.pT,clx.Tv,NMC);
    
    
    % main loop 
    [dpc,sys] = dpc_ini(sys,dat,opt.T,sys.rho_star);
    clx.Ac = [eye(dpc.N) -eye(dpc.N); -eye(dpc.N) -eye(dpc.N)];
    clx.bc = zeros(2*dpc.N,1);
    clx.TN = opt.T/dpc.N;
    clx.nsvars = nsvars;
    
    
    for j = 1:NMC
        [clx,sys,dpc,exe_time,...
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
        JFCE,JFCE_u,JFCE_y,u_FCE,y_FCE,yFCE_TR,yFCE_KF,beta_FCE,...
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
        J23_bar,Ju23_bar,Jy23_bar,u23_bar,y23_bar,y23_bar_TR,...
        y23_bar_KF,n2_gamma12_23_bar,gamma3_23_bar,...
        J1_bar,Ju1_bar,Jy1_bar,u1_bar,y1_bar,y1_bar_TR,y1_bar_KF,...
        JFCE,JFCE_u,JFCE_y,u_FCE,y_FCE,yFCE_TR,yFCE_KF,beta_FCE,...
        Jthm3_tuned,Ju_thm3_tuned,Jy_thm3_tuned,u_thm3_tuned,...
        y_thm3_tuned,y_thm3_tuned_TR,...
        y_thm3_tuned_KF,n2_gamma12_thm3_tuned,beta_thm3_tuned,...
        J2_tuned,Ju2_tuned,Jy2_tuned,u2_tuned,y2_tuned,y2_tuned_TR,...
        y2_tuned_KF,n2_gamma12_2_tuned,beta2_tuned,SD2abs_actual,...
        J3_tuned,Ju3_tuned,Jy3_tuned,u3_tuned,y3_tuned,y3_tuned_TR,...
        y3_tuned_KF,n2_gamma12_3_tuned,gamma3_tuned,beta3_tuned,...
        SD3abs_actual,dat,opt,gra,rho,idx1,idx2,idx_1,idx_2,...
        rho_search_j,Lpts,NMC_inner,j);
    
        exe_times{j} = exe_time;
        if savefig_flag && mod(j,10) == 5
            save(['j' num2str(j) '_model' num2str(model)...
                '_SNR' num2str(snr_target_dB) '_wrkspc.mat'])
        end
    end
    
    
    %% MC averaging
    fprintf('\n__________________________\n')
    fprintf('MC averaging...\n\n')
    
    % beta bar found and related costs
    J1_ave = J1_ave / NMC;
    J2_ave = J2_ave / NMC;
    J3_ave = J3_ave / NMC;
    J23_ave = J23_ave / NMC;
    [J2_min,i2_min_bar] = min(J2_ave);
    [J3_min,i3_min_bar] = min(J3_ave);
    J23_min = min(min(J23_ave));
    [i22_min_bar,i33_min_bar] = find(J23_ave == J23_min);
    J1_min = min(min(J1_ave));
    [il1_min_bar,il2_min_bar] = find(J1_ave == J1_min);
    best_beta2_bar = clx.beta2_(i2_min_bar);
    best_beta3_bar = clx.beta3_(i3_min_bar);
    best_beta22_bar = clx.beta2_(i22_min_bar(1));
    best_beta33_bar = clx.beta3_(i33_min_bar(1));
    best_lmb1_bar = clx.lmb1_(il1_min_bar(1));
    best_lmb2_bar = clx.lmb2_(il2_min_bar(1));
    MC_beta2_bar = median(beta2_bar);
    MC_beta3_bar = median(beta3_bar);
    MC_beta22_bar = median(best_beta22_bar);
    MC_beta33_bar = median(best_beta33_bar);
    MC_best_lmb1_bar = median(best_lmb1_bar);
    MC_best_lmb2_bar = median(best_lmb2_bar);
     
    fprintf('\n**************************\n')
    fprintf('Found betas on grid:\n')
    fprintf(['best_beta2_bar = ' num2str(best_beta2_bar) '\n'])
    fprintf(['best_beta3_bar = ' num2str(best_beta3_bar) '\n'])
    fprintf(['best_beta22_bar = ' num2str(best_beta22_bar) '\n'])
    fprintf(['best_beta33_bar = ' num2str(best_beta33_bar) '\n'])
    fprintf(['MC_beta2_bar = ' num2str(MC_beta2_bar) '\n'])
    fprintf(['MC_beta3_bar = ' num2str(MC_beta3_bar) '\n'])
    fprintf(['MC_beta22_bar = ' num2str(MC_beta22_bar) '\n'])
    fprintf(['MC_beta33_bar = ' num2str(MC_beta33_bar) '\n'])
    fprintf(['MC_best_lmb1_bar = ' num2str(MC_best_lmb1_bar) '\n'])
    fprintf(['MC_best_lmb2_bar = ' num2str(MC_best_lmb2_bar) '\n'])
    fprintf('**************************\n')
    
    % average execution times
    ave_exe_times = zeros(18,1);
    for j = 1:NMC
        ave_exe_times(1) = ave_exe_times(1) + exe_times{j}.KF;
        ave_exe_times(2) = ave_exe_times(2) + exe_times{j}.b2and3g;
        ave_exe_times(3) = ave_exe_times(3) + exe_times{j}.b23g;
        ave_exe_times(4) = ave_exe_times(4) + exe_times{j}.l12g;
        ave_exe_times(5) = ave_exe_times(5) + exe_times{j}.b2;
        ave_exe_times(6) = ave_exe_times(6) + exe_times{j}.b3;
        ave_exe_times(7) = ave_exe_times(7) + exe_times{j}.b23;
        ave_exe_times(8) = ave_exe_times(8) + exe_times{j}.l12;
        ave_exe_times(9) = ave_exe_times(9) + exe_times{j}.O;
        ave_exe_times(10) = ave_exe_times(10) + exe_times{j}.S;
        ave_exe_times(11) = ave_exe_times(11) + exe_times{j}.b2t;
        ave_exe_times(12) = ave_exe_times(12) + exe_times{j}.b3t;
        ave_exe_times(13) = ave_exe_times(13) + exe_times{j}.trn_others;
        ave_exe_times(14) = ave_exe_times(14) + exe_times{j}.trn_PBSIDopt;
        ave_exe_times(15) = ave_exe_times(15) + exe_times{j}.trn_S;
        ave_exe_times(16) = ave_exe_times(16) + exe_times{j}.trn_23online;
        ave_exe_times(17) = ave_exe_times(17) + exe_times{j}.trn_DeePC;
    end
    ave_exe_times(18) = trn_rho;
    ave_exe_times = ave_exe_times / NMC;
    table.training = [-1 ave_exe_times(13)+ave_exe_times(18) ...
        ave_exe_times(13)+ave_exe_times(18) ...
        ave_exe_times(13)+ave_exe_times(18) ...
        ave_exe_times(17)+ave_exe_times(18) ...
        ave_exe_times(14)+ave_exe_times(18) ...
        ave_exe_times(15)+ave_exe_times(18) ...
        ave_exe_times(16)+ave_exe_times(18) ...
        ave_exe_times(16)+ave_exe_times(18)];
    table.offlsearch = [-1 ave_exe_times(2)/2 ave_exe_times(2)/2 ...
        ave_exe_times(3) ave_exe_times(4) -1 -1 -1 -1];
    table.optimization = [ave_exe_times(1) ave_exe_times(5) ...
        ave_exe_times(6) ave_exe_times(7) ave_exe_times(8) ...
        ave_exe_times(9) ave_exe_times(10) ave_exe_times(11) ...
        ave_exe_times(12)];
    fprintf(['Time MPC: ' num2str(ave_exe_times(1)) ' s \n'])
    fprintf(['Time b2and3g: ' num2str(ave_exe_times(2)) ' s \n'])
    fprintf(['Time b23g: ' num2str(ave_exe_times(3)/Lpts_r) ' s \n'])
    fprintf(['Time l12g: ' num2str(ave_exe_times(4)/Lpts_r_lmb) ' s \n'])
    fprintf(['Time b2: ' num2str(ave_exe_times(5)) ' s \n'])
    fprintf(['Time b3: ' num2str(ave_exe_times(6)) ' s \n'])
    fprintf(['Time b23: ' num2str(ave_exe_times(7)) ' s \n'])
    fprintf(['Time l12: ' num2str(ave_exe_times(8)) ' s \n'])
    fprintf(['Time FCE: ' num2str(ave_exe_times(9)) ' s \n'])
    fprintf(['Time thm3: ' num2str(ave_exe_times(10)) ' s \n'])
    fprintf(['Time b2t: ' num2str(ave_exe_times(11)) ' s \n'])
    fprintf(['Time b3t: ' num2str(ave_exe_times(12)) ' s \n'])

    MC_J_MPC = mean(J_MPC,Ls(J_MPC));
    MC_J2_grid = mean(J2_grid,Ls(J2_grid));
    MC_J3_grid = mean(J3_grid,Ls(J3_grid));
    MC_J2_bar = mean(J2_bar);
    MC_J3_bar = mean(J3_bar);
    MC_J1_bar = mean(J1_bar,Ls(J1_bar)); 
    MC_J23_bar = mean(J23_bar,Ls(J23_bar));        
    
    % tuned stuff  
    MC_Jthm3_tuned = mean(Jthm3_tuned,Ls(Jthm3_tuned));            % scalar
    MC_JFCE = mean(JFCE,Ls(JFCE));                                 % scalar
    MC_J2_tuned = mean(J2_tuned,Ls(J2_tuned));                     % scalar
    MC_J3_tuned = mean(J3_tuned,Ls(J3_tuned));                     % scalar
    MC_SD2_actual = mean(SD2abs_actual,Ls(SD2abs_actual));         % Tv
    mean_MC_SD2_actual = mean(MC_SD2_actual,Ls(MC_SD2_actual));    % scalar
    MC_beta2_tuned = mean(beta2_tuned,Ls(beta2_tuned));            % Tv
    MC_beta_thm3 = mean(beta_thm3_tuned,Ls(beta_thm3_tuned));      % Tv
    MC_beta_FCE = mean(beta_FCE,Ls(beta_FCE));                     % Tv
    mean_MC_beta2_tuned = mean(MC_beta2_tuned,Ls(MC_beta2_tuned)); % scalar
    MC_beta3_tuned = mean(beta3_tuned,Ls(beta3_tuned));            % Tv
    mean_MC_beta3_tuned = mean(MC_beta3_tuned,Ls(MC_beta3_tuned)); % scalar
    MC_SD3_actual = mean(SD3abs_actual,Ls(SD3abs_actual));         % Tv
    mean_MC_SD3_actual = mean(MC_SD3_actual,Ls(MC_SD3_actual));    % scalar
    
    
    % error('Checkpoint: MC averaging done')
    

    %% showing final numerical results
    fprintf('\n**************************\n')
    fprintf('Final results:\n')
    fprintf('________\n')
    fprintf(['J_GT = ' num2str(sol_GT.J) '\n'])
    fprintf(['MC_J_MPC = ' num2str(MC_J_MPC) '\n'])
    fprintf(['MC J2 bar = ' num2str(MC_J2_bar) '\n'])
    fprintf(['MC J3 bar = ' num2str(MC_J3_bar) '\n'])
    fprintf(['MC J23 bar = ' num2str(MC_J23_bar) '\n'])
    fprintf(['MC JDeePC bar = ' num2str(MC_J1_bar) '\n'])
    fprintf(['MC JFCE tuned = ' num2str(MC_JFCE) '\n'])
    fprintf(['MC Jthm3 tuned = ' num2str(MC_Jthm3_tuned) '\n'])
    fprintf(['MC J2 tuned = ' num2str(MC_J2_tuned) '\n'])
    fprintf(['MC J3 tuned = ' num2str(MC_J3_tuned) '\n'])
    fprintf('________\n')
    % display(mean_MC_SD2_actual)
    % display(mean_MC_SD3_actual)
    % display(mean_MC_beta2_tuned)
    % display(mean_MC_beta3_tuned)
    fprintf('**************************\n')
    
    
    % error('Checkpoint: MC final results')
    

    %% Producing plots
    if savefig_flag
    fprintf('\n____________________________\n')
    fprintf('Producing plots ...\n')
    fprintf('____________________________\n')
    
    
    
    xempty = [];
    parfor j = 1:NMC
        xempty = strcat(xempty,'-');
    end
    xempty = xempty';

    % my_boxplot(x,y,constants,legs,xlab,ylab,titl,Ylog)

    % All final CL performances
    gra.ftsz = 28;
    xcell = {'$J_{MPC}$', '$\bar{J}_{\gamma_{2}}$',... 
        '$\bar{J}_{\gamma_{3}}$',... 
        '$\bar{J}_{\gamma_{23}}$','$\bar{J}_{DeePC}$',... %
          '$\hat{J}_{FCE}$', '$\hat{J}_{thm3}$',... 
         '$\hat{J}_{\gamma_{2}}$', '$\hat{J}_{\gamma_{3}}$'};
    gra.pos = [100 100 800 700];
    
    my_boxplot(gra,xcell,[J_MPC(1:NMC); 
        J2_bar(1:NMC); J3_bar(1:NMC); J23_bar(1:NMC); J1_bar(1:NMC); ...
        JFCE(1:NMC); Jthm3_tuned(1:NMC); J2_tuned(1:NMC); ...
        J3_tuned(1:NMC)],...
        [],...
        {},...
        'northeast',...
        '','','Closed-loop performance',0,0);
    if SNRs_dB == 6
        ylim([0 5]);
        yticks(0:5);
        if model == 3
            number_of_instabilities = sum(J3_tuned > 5)
        end
    elseif SNRs_dB == 20
        ylim([0 0.25]);
        yticks(0:0.05:0.25);
    end
    
    
    if 1
    % All final CL control costs
    xcell = {'$J_{u,MPC}$', '$\bar{J}_{u,2}$', '$\bar{J}_{u,3}$',...
        '$\bar{J}_{u,23}$', '$\bar{J}_{u,DeePC}$', '$\hat{J}_{u,FCE}$',...
        '$\hat{J}_{u,thm3}$', '$\hat{J}_{u,2}$', '$\hat{J}_{u,3}$', };
    my_boxplot(gra,xcell,[J_MPC_u; Ju2_bar; Ju3_bar; Ju23_bar; Ju1_bar;...
        JFCE_u; Ju_thm3_tuned; Ju2_tuned; Ju3_tuned],...
        [],{},'','','',...
        'Closed-loop control effort',0,0);
    
    % All final CL relative tracking errors (normalized by yr)
    xcell = {'$J_{y,MPC}$', '$\bar{J}_{y,2}$', '$\bar{J}_{y,3}$',...
        '$\bar{J}_{y,23}$', '$\bar{J}_{y,DeePC}$', '$\hat{J}_{y,FCE}$',...
        '$\hat{J}_{y,thm3}$', '$\hat{J}_{y,2}$', '$\hat{J}_{y,3}$'};
    my_boxplot(gra,xcell,[J_MPC_y; Jy2_bar; Jy3_bar; Jy23_bar; Jy1_bar;...
        JFCE_y; Jy_thm3_tuned; Jy2_tuned; Jy3_tuned],...
        [],{},'','','',...
        'Closed-loop relative tracking error (normalized by $y_r$)',0,0);
    
    % Distributions of beta_bar
    xcell = {'$\bar{\beta}_{2}$', '$\bar{\beta}_{3}$'};
    my_boxplot(gra,xcell,[beta2_bar; beta3_bar],...
        [],{},'','','',...
        'Distributions of $\bar{\beta}_{2,3}$',1,0);
    
    % beta2 and beta3 comparison for minimum CL average costs
    if ~sys.kind_of_setup
        gra.ftsz = 20;
        gra.pos = [0 100 1100 500];
        figure('position',gra.pos)
        fig = subplot(2,1,1);
        xlab = '$\beta_2$';
        ylab = '$J_{AV,2}(\beta_2)$';
        titl = 'Average closed-loop performances vs $\beta_2$';
        ave_cost_fig(gra,fig,clx.beta2_,best_beta2_bar,J2_ave,J2_min,...
            xlab,ylab,titl,[1 0 0]);
        fig = subplot(2,1,2);
        xlab = '$\beta_3$';
        ylab = '$J_{AV,3}(\beta_3)$';
        titl = 'Average closed-loop performances vs $\beta_3$';
        ave_cost_fig(gra,fig,clx.beta3_,best_beta3_bar,J3_ave,J3_min,...
            xlab,ylab,titl,[0 0 1]);
    end
    
    
    
    %% tracking errors
    if ~sys.kind_of_setup
        gra.ftsz = 30;
        gra.pos = [100 100 600 600];
        my_col = [128,0,128]/255;
        pp = 1:sys.p;
        track_err_fig1(gra,clx,NMC,y0(pp,:,1:NMC), ...
            '$t$','$y_{MPC}(t)$','',my_col)
        track_err_fig1(gra,clx,NMC,y1_bar_TR(pp,:,1:NMC), ...
            '$t$','$y_{DeePC}(t)$','',my_col)
        my_col = [1 0 0];
        track_err_fig1(gra,clx,NMC,y2_bar_TR(pp,:,1:NMC), ...
            '$t$','$\bar{y}_2(t)$','',my_col)
        track_err_fig1(gra,clx,NMC,y2_tuned_TR(pp,:,1:NMC), ...
            '$t$','$\hat{y}_2(t)$','',my_col)
        my_col = [0 0 1];
        track_err_fig1(gra,clx,NMC,y3_bar_TR(pp,:,1:NMC), ...
            '$t$','$\bar{y}_3(t)$','',my_col)
        track_err_fig1(gra,clx,NMC,y3_tuned_TR(pp,:,1:NMC), ...
            '$t$','$\hat{y}_3(t)$','',my_col)
        my_col = [1 0 1];
        track_err_fig1(gra,clx,NMC,y23_bar_TR(pp,:,1:NMC), ...
            '$t$','$\bar{y}_{23}(t)$','',my_col)
        my_col = [0 1 0];
        track_err_fig1(gra,clx,NMC,yFCE_TR(pp,:,1:NMC), ...
            '$t$','$\hat{y}_{FCE}(t)$','',my_col)
        my_col = [1 1 0];
        track_err_fig1(gra,clx,NMC,y_thm3_tuned_TR(pp,:,1:NMC), ...
            '$t$','$\hat{y}_{thm3}(t)$','',my_col)
    else
        %% tracking errors
        % NOTE
        % use set(gcf,'renderer','Painters')
        % to save high quality figures
        gra.ftsz = 30;
        gra.pos = [100 100 600 600];
        my_col = [128,0,128]/255;
        pp = 1:sys.p;
        % track_err_fig1(gra,clx,NMC,y1(pp,:,1:NMC),...
        % '$t$','$y_{MPC}(t)$','',my_col)
        my_col = [1 0 0];
        perc = 50;
        gra.j_star = getPercIndex(Jy1_bar,perc);
        % DEACTIVATE THE ERROR BELOW 
        error('Make sure to configure track_err_fig3.m as required')
        track_err_fig3(gra,clx,NMC,y1_bar_TR(pp,:,1:NMC),...
            '$t$','$y_{DeePC}(t)$','',my_col)
        my_col = [1 0 0];
        % track_err_fig1(gra,clx,NMC,y2_bar_TR(pp,:,1:NMC),...
        % '$t$','$\bar{y}_2(t)$','',my_col)
        % track_err_fig1(gra,clx,NMC,y2_tuned_TR(pp,:,1:NMC),...
        % '$t$','$\hat{y}_2(t)$','',my_col)
        my_col = [0 0 1];
        % track_err_fig1(gra,clx,NMC,y3_bar_TR(pp,:,1:NMC),...
        % '$t$','$\bar{y}_3(t)$','',my_col)
        % track_err_fig1(gra,clx,NMC,y3_tuned_TR(pp,:,1:NMC),...
        % '$t$','$\hat{y}_3(t)$','',my_col)
        my_col = [1 0 1];
        %track_err_fig1(gra,clx,NMC,y23_bar_TR(pp,:,1:NMC),'$t$',...
        % '$\bar{y}_{23}(t)$','',my_col)
        my_col = [0 0 1];
        gra.j_star = getPercIndex(Jyo,perc);
        track_err_fig3(gra,clx,NMC,yo_TR(pp,:,1:NMC),...
            '$t$','$y_{FCE}(t)$','',my_col)
        my_col = [1 1 0];
        % track_err_fig1(gra,clx,NMC,y_thm3_tuned_TR(pp,:,1:NMC),...
        % '$t$','$\hat{y}_{S}(t)$','',my_col)
        %%
        if 0
            pp = 1:sys.p;
            my_col = [1 0 0];
            weird_trajectories(gra,clx,NMC,y1_bar_TR(pp,:,1:NMC),...
                '$t$','$y_{DeePC}(t)$','',my_col)
        end
    end
    
    %%
    mm = 1:sys.m;
    track_err_fig1(gra,clx,NMC,u0(mm,:,:), ...
        '$t$','$u_{MPC}(t)$','',my_col)
    my_col = [1 0 0];
    track_err_fig1(gra,clx,NMC,u2_bar(mm,:,:), ...
        '$t$','$\bar{u}_2(t)$','',my_col)
    track_err_fig1(gra,clx,NMC,u2_tuned(mm,:,:), ...
        '$t$','$\hat{u}_2(t)$','',my_col)
    my_col = [0 0 1];
    track_err_fig1(gra,clx,NMC,u3_bar(mm,:,:), ...
        '$t$','$\bar{u}_3(t)$','',my_col)
    track_err_fig1(gra,clx,NMC,u3_tuned(mm,:,:), ...
        '$t$','$\hat{u}_3(t)$','',my_col)
    my_col = [1 0 1];
    track_err_fig1(gra,clx,NMC,u23_bar(mm,:,:), ...
        '$t$','$\bar{u}_{23}(t)$','',my_col)
    my_col = [0 1 0];
    track_err_fig1(gra,clx,NMC,u_FCE(mm,:,:), ...
        '$t$','$\hat{u}_{FCE}(t)$','',my_col)
    my_col = [1 1 0];
    track_err_fig1(gra,clx,NMC,u_thm3_tuned(mm,:,:), ...
        '$t$','$\hat{u}_{thm3}(t)$','',my_col)
    %%
    
    end

    end
    
end


if savefig_flag
    fig_name = ...
        ['S' num2str(sys.kind_of_setup) '_model' num2str(model) '_SNR' ...
        num2str(snr_target_dB) '.fig'];
    savefig(fig_counter:fig_counter+nfigs-1,fig_name);
    save(['S' num2str(sys.kind_of_setup) '_model' num2str(model) '_SNR' ...
        num2str(snr_target_dB) '_wrkspc.mat'])
end

end
end

%% End of program
fprintf('DONE\n')
main_time = toc(main_time);
cvx_sol_share = cvx_sol_time/main_time;
cvx_sol_self_time = cvx_sol_time/cvx_sol_counter;
main_time
cvx_sol_time
cvx_sol_share
cvx_sol_counter
cvx_sol_self_time


function j_star = getPercIndex(J,perc)
    perc = perc/100;
    [~,idx] = sort(J);
    j_star = idx(ceil(perc*length(J)));
end
