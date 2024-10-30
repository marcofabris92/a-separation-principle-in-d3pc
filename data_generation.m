% This function is used for two purposes:
% - adjust some system parameters, such as those related to the SNR ratios
%   and covariance matrices, and compute some other parameters, such as
%   those related to Kalman filters used in the computation of the MPC 
%   solution; this functionality completes what could not be done yet in 
%   model_selection.m (see part 1) of this piece of code)
% - generate synthetic data (stored in dat), i.e. the training data sets
%   (see part 2) and 3) of this piece of code); here the struct dat is 
%   initialized (exept for the field Ndata)
% Given
% - the underlying LTI system parameters contained in sys
% - the number of data to be contained in the data sets being created,
%   stored in dat.Ndata
% - the optimization parameters and references contained in opt
% - the standard deviation sigma_u for the input training data
% - the desired SNR, stored in SNR_dB
% - the flag clr, if it is set to 1 then filtering is applied on the input
%   u_long, which is use to generate all the data sets
% - the total number of Monte Carlo runs NMC

% Invoked by:
% - MAIN.m, to initialize the synthetic data sets for the simulations
% Invokes:
% - dynamics.m, to construct input/output sequences
% - filter_design.m, whenever filtering on u_long is desired
% - sample_var(), to compute the sample variance of a sequence of data
%   (see below)


function [sys,dat] = data_generation(sys,dat,opt,sigma_u,SNR_dB,clr,NMC)

x0 = zeros(sys.n,1);
u_long = randn(sys.m,max(10*dat.Ndata*sys.p,1e6));   

%% 1) adjusting sys.Q, sys.R, sys.S; computing Kalman parameters 
if sys.flag_noise
    fprintf('Adjusting SNR ...\n\n')
    
    % getting the sqrtm of noise covariance matrix
    if isempty(sys.K)
        rk_V = rank([sys.Q sys.S; sys.S' sys.R]);
        [U_V12,S_V12,~] = svd([sys.Q sys.S; sys.S' sys.R]);
        sys.V12 = U_V12(:,1:rk_V)*sqrtm(S_V12(1:rk_V,1:rk_V));
    else
        sys.V12 = [sys.K; eye(sys.p)]*sqrtm(sys.R);
    end

    % adjustment of SNR ratios
    y_stoch_comp_long = dynamics(sys,x0,zeros(size(u_long)),[],[]);
    sys.flag_noise = 0;
    y_det_long = dynamics(sys,x0,u_long,[],[]);
    sys.flag_noise = 1;
    snr_0 = 1+sample_var(y_det_long)/sample_var(y_stoch_comp_long);
    snr_target = 10^(SNR_dB/10);
    snr_ratio = (snr_0-1)/(snr_target-1);

    % rescaling state and output noise covariances
    sys.Q = snr_ratio * sys.Q; % adjusted Q
    sys.R = snr_ratio * sys.R; % adjusted R
    sys.S = snr_ratio * sys.S; % adjusted S

    % computing the Kalman update matrices
    sys.Z = sys.S*(sys.R^-1);
    sys.F = sys.A-sys.Z*sys.C;
    sys.B_bar = sys.B-sys.Z*sys.D;
    sys.Q_tilde = sys.Q-sys.Z*sys.S';
    if ~isempty(sys.K)
        sys.Q_tilde = zeros(sys.n);
    end
    
    % getting the sqrtm of noise covariance matrix
    if isempty(sys.K)
        rk_V = rank([sys.Q sys.S; sys.S' sys.R]);
        [U_V12,S_V12,~] = svd([sys.Q sys.S; sys.S' sys.R]);
        sys.V12 = U_V12(:,1:rk_V)*sqrtm(S_V12(1:rk_V,1:rk_V));
    else
        sys.V12 = [sys.K; eye(sys.p)]*sqrtm(sys.R);
        V = sys.V12*sys.V12';
        sys.Q = V(1:sys.n,1:sys.n);
    end
    
    % getting Kalman stationary matrices (gain and innovation variance)
    [sys.P_inf,~,sys.K_inf] = dare(sys.F',sys.C',sys.Q_tilde,sys.R);
    sys.K_inf = sys.K_inf';
    sys.Lambda_inf = sys.C*sys.P_inf*sys.C'+sys.R;
    
    % getting innovation form matrices
    if isempty(sys.K)
        sys.K = sys.K_inf+sys.Z;
        sys.Lambda = sys.Lambda_inf;
    else
        sys.Lambda = sys.R;
    end
    
    % compute the variance of the stochastic part of the total cost
    Lam12 = sqrtm(sys.Lambda);
    sys.stoch_cost_Var = trace(Lam12'*opt.Q*Lam12);
    
    % computation of Hs
    sys.Hs = eye(opt.pT);
    for low_diag = 0:opt.T-2
        ii = (1:sys.p)+low_diag*sys.p;
        new_entry = sys.Gamma(ii,:)*sys.K;
        for col = low_diag:opt.T-2
            ii = (1:sys.p)+(col+1)*sys.p;
            jj = (1:sys.p)+(col-low_diag)*sys.p;
            sys.Hs(ii,jj) = new_entry;
        end
    end
    
    % computation of inverse and normalized Hs
    sys.HsLam_inv = pinv(sys.Hs*kron(opt.I_T,sqrtm(sys.Lambda)));
    
    % computation of auxiliary matrices used in ol_MPC.m
    Upsilon = (opt.TR+sys.Hd'*opt.TQ*sys.Hd)^-1;
    sys.Upsilon_u = Upsilon*opt.TR;
    sys.Upsilon_y = Upsilon*sys.Hd'*opt.TQ;
    
    % assigning state of the predictor and prediction error variance
    % default choice: it will be arbitrarily overridden in dpc_ini.m!
    sys.P = dlyap(sys.F,sys.Q);
    sys.x_ini__ini_1 = zeros(sys.n,1);
    sys.P_ini__ini_1 = sys.P;
    
    % choosing the initial state to generate data
    rk_P = rank(sys.P);
    [U_P12,S_P12,~] = svd(sys.P);
    sys.P12 = U_P12(:,1:rk_P)*sqrtm(S_P12(1:rk_P,1:rk_P));
    x0 = 0*sys.P12*randn(size(sys.P12,2),1);
    
    % choosing the initial state of the underlying problem
    sys.x0 = 0*sys.P12*randn(size(sys.P12,2),1);
    

end
    
%% 2) initializing dat.u and dat.y (noiseless data set)
if clr
    HLP = filter_design(1,1/pi);
    my_filter = tf(HLP.Numerator,[1 zeros(1,length(HLP.Numerator)-1)],...
        sys.Ts);
    l_tr = round(1.2*length(HLP.Numerator));
    dat.u = lsim(my_filter,u_long(:,1:dat.Ndata+l_tr))';
    dat.u = dat.u(:,end-dat.Ndata+1:end);
else
    dat.u = u_long(:,1:dat.Ndata);
end
flag_noise_copy = sys.flag_noise;
sys.flag_noise = 0;
dat.y = dynamics(sys,x0,dat.u,[],[]);
sys.flag_noise = flag_noise_copy;

%% 3) getting training dataset (NMC noisy data sets)
dat.u_noisy = zeros(sys.m,dat.Ndata,NMC);
dat.y_noisy = zeros(sys.p,dat.Ndata,NMC);
for j = 1:NMC
    dat.u_noisy(:,:,j) = sigma_u*randn(sys.m,dat.Ndata);
    dat.y_noisy(:,:,j) = dynamics(sys,x0,dat.u_noisy(:,:,j),[],[]);
end

end



% This function is used to compute the sample matrix variance of a 
% sequence x, which is returned as V_hat
% Given
% - the sequence x

% Invokes: none
% Invoked by:
% - data_generation.m, to adjust the SNR as desired


function [V_hat] = sample_var(x)

% rows of x: features
% columns of x: observations (observations >> features is assumed)

V_hat = 0;
mu = mean(x,2);
T_long = size(x,2);
parfor t = 1:T_long
    V_hat = V_hat + norm(x(:,t)-mu)^2;
end
V_hat = V_hat / (T_long-size(x,1));

end