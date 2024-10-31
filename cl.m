% This function starts the closed-loop to run the control being tested on
% the selected dynamics. It proceeds by computing the solution, step by
% step until it reaches step Tv. It returns the struct sol containing all
% the results related to the solution corresponding to a specific 
% predictive control approach.
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given data set contained in dpc
% - type: an integer that determines which appraoch should be used (see
%   below), at the beginning of cl()
% - MPC_comp, a flag = 1 if the comparison with standard MPC is required, 
%   given the same input signal computed through any other data-driven
%   approach
% - beta, containing the value of beta for DDPC approaches or the values of
%   lambda for DeePC (all these are just already selected scaling factors 
% for the regularizers); beta = [] if online tuning is required

% Invoked by:
% - get_GT.m, to launch the closed-loop dynamics for the computation of the
%   noise-free solution
% - innerMC_grids.m, to perform offline grid optimizations
% - outerMC.m, to execute the closed-loop dynamics and compute the
%   solutions of all the predictive control approaches being compared
% Invokes:
% - dynamics.m, to push forward each step of the closed-loop
% - ol(), to launch all the different predictive control approaches
%   (see below)
% - cl_cost_u.m, to retrieve information of the control cost
% - cl_cost_y.m, to retrieve information of the tracking cost
% - cl_cost.m, to compute the performance index of interest based on the
%   weighted sum of control and tracking costs


function [sol] = cl(clx,dpc,type,MPC_comp,beta)

% type
% 0: noise-free ground truth (used just for sanity check)
% 1: MPC solution based on the knowledge of the "true" model
% 2: gamma-DDPC with gamma_2 only in the regularization term (see [9])
% 3: gamma-DDPC with gamma_3 only in the regularization term (see [9])
% 4: gamma-DDPC with both gamma_2 and gamma_3 in the regularization term
%    (see [7])
% 5: DeePC in (25) of [15]
% 6: FCE (proposed method in this paper)
% 7: thm3 (proposed method in this paper "Theorem 3", less accurate)


%% COMMON INITIALIZATION
% unpacking variables
sys = clx.sys;
opt = clx.opt;
Tv = clx.Tv;
ddpc_approach = (type == 0 || type > 1);
KF_approach = (type == 1);
mm = 1:sys.m;
pp = 1:sys.p;

% solution initialization
sol.uf = zeros(opt.mT,Tv);
sol.yf_hat = zeros(opt.pT,Tv);
sol.yf = zeros(opt.pT,Tv);

%% SPECIFIC INITIALIZATION
% initialization of the past trajectories
xini = dpc.xini;                       
vx_true = [];
vy_true = [];
if ddpc_approach
    uini = dpc.vec_uini;
    yini = dpc.vec_yini;
    if type ~= 6
        sol.n2gamma12 = zeros(1,Tv);
    end
end
if KF_approach
    prd.x_t__t_1 = sys.x_ini__ini_1;    % x0|-1
    prd.P_t__t_1 = sys.P_ini__ini_1;    % P0|-1
end
if ddpc_approach && MPC_comp
    sol.yf_hat_star = zeros(opt.pT,Tv);
    prd_KF_comp.x_t__t_1 = sys.x_ini__ini_1;    % x0|-1
    prd_KF_comp.P_t__t_1 = sys.P_ini__ini_1;    % P0|-1
end

% initializing interesting quantities depending on "type"
eb = isempty(beta);
switch type
    case 0
        sys.flag_noise = 0;
    case 1
    case 2
        if eb
            sol.beta2_tuned = zeros(1,Tv);
            sol.actual_SD2 = zeros(1,Tv);
        end
    case 3
        sol.gamma3 = zeros(opt.pT,Tv);
        if eb
            sol.beta3_tuned = zeros(1,Tv);
            sol.actual_SD3 = zeros(1,Tv);
        end
    case 4
        sol.gamma3 = zeros(opt.pT,Tv);
    case 5
        prd.lmb1 = beta(1);
        prd.lmb2 = beta(2);
        prd.options = optimoptions('fmincon');
        prd.options.SpecifyObjectiveGradient = true;
        prd.options.ConstraintTolerance = 1e-4;
        prd.options.MaxFunctionEvaluations = 2e1;
        if prd.lmb1 < +Inf
            HH = [...
            (dpc.Yf'*opt.TQ*dpc.Yf+dpc.Uf'*opt.TR*dpc.Uf + ...
            prd.lmb1*dpc.IPi2)                      zeros(dpc.N);
            zeros(dpc.N)                            zeros(dpc.N)];
            prd.options.HessianFcn = @(x,lambda) HH;
        end
        prd.options.Display = 'off';
        L22_inv = dpc.L22^-1;
        L33_inv = dpc.L33^-1;
    case 6 
        if eb
            sol.beta_FCE = zeros(1,Tv);
        end
    case 7
        if eb
            sol.beta_thm3 = zeros(1,Tv);
        end
end

if type == 6
    sol.J_yu = 0;
    sol.J_reg = 0;
end

% main closed-loop
for t = 1:Tv
    %% SOLVING
    prd.ur = opt.vec_ur((t-1)*sys.m+1 : (t+opt.T-1)*sys.m);
    prd.yr = opt.vec_yr((t-1)*sys.p+1 : (t+opt.T-1)*sys.p);
    if ddpc_approach
        if type == 6
            prd.zini = zeros(size(uini,1)+size(yini,1),1);
            for i = 0:dpc.Tini-1 
                ii = (sys.m+sys.p)*i+1 : (sys.m+sys.p)*(i+1);
                ip = sys.p*i+1 : sys.p*(i+1);
                im = sys.m*i+1 : sys.m*(i+1);
                prd.zini(ii) = [yini(ip); uini(im)];
            end
        else
            prd.zini = [uini(end-dpc.Tini*sys.m+1:end);...
                yini(end-dpc.Tini*sys.p+1:end)];
            prd.gamma1 = dpc.pinv_L11*prd.zini;
        end
        if type == 4
            prd.c3 = opt.TQ*(prd.yr-dpc.L31*prd.gamma1);
            prd.c2 = dpc.L32'*prd.c3+...
                dpc.L22'*opt.TR*(prd.ur-dpc.L21*prd.gamma1);
            prd.c3 = dpc.L33'*prd.c3;
        end
        if type == 5
            gamma2 = L22_inv*(prd.ur-dpc.L21*prd.gamma1);
            gamma3 = L33_inv*(prd.yr-dpc.L31*prd.gamma1-dpc.L32*gamma2);
            prd.g0 = dpc.Q_inv*[prd.gamma1; gamma2; gamma3];
        end
    end
    prd = ol(clx,dpc,xini,type,prd,beta,t);
    if sys.flag_noise
        vx_true = sys.vx(:,t:t+opt.T-1);
        vy_true = sys.vy(:,t:t+opt.T-1);
    end
    u_true = reshape(prd.u,[sys.m,opt.T]);
    [y_true,x_true] = dynamics(sys,xini,u_true,vx_true,vy_true);
    %% COMMON ASSIGNMENT
    % assigning solutions
    sol.uf(:,t) = prd.u;
    sol.yf_hat(:,t) = prd.y;
    sol.yf(:,t) = reshape(y_true,[opt.pT,1]);
    
    %% SPECIFIC ASSIGNMENT
    % assigning interesting quantities depending on "type"
    if ddpc_approach
        if type ~= 5 && type ~= 6 && type ~= 7
            sol.n2gamma12(t) = norm([prd.gamma1; prd.gamma2])^2;
        end
        if MPC_comp
            prd_KF_comp.ur = prd.ur;
            prd_KF_comp.yr = prd.yr;
            prd_KF_comp = ol_MPC(clx,xini,prd_KF_comp,prd.u,t);
            sol.yf_hat_star(:,t) = prd_KF_comp.y;
        end
    end
    switch type
        case 0
        case 1
        case 2
            if eb
                sol.beta2_tuned(t) = prd.beta;
                sol.actual_SD2(t) = SD_2(clx,prd);
            end
        case 3
            sol.gamma3(:,t) = prd.gamma3;
            if eb
                sol.beta3_tuned(t) = prd.beta;
                sol.actual_SD3(t) = SD_3(clx,prd);
            end
        case 4
            sol.gamma3(:,t) = prd.gamma3;
        case 6
            if eb
                sol.beta_FCE(t) = prd.beta;
            end
        case 7
            if eb
                sol.beta_thm3(t) = prd.beta;
            end
    end
    %% UPDATE
    xini = x_true(:,2);
    if ddpc_approach
        uini = [uini(sys.m+1:end); prd.u(mm)];
        yini = [yini(sys.p+1:end); y_true(:,1)];
    end
    if type == 6
        sol.J_yu = sol.J_yu + prd.yu_cost;
        sol.J_reg = sol.J_reg + prd.reg_cost;
    end
end

% computation of the closed-loop performance
sol.Ju = cl_cost_u(opt,sol.uf(mm,:))/Tv;
sol.Jy = cl_cost_y(opt,sol.yf(pp,:))/opt.yrTvnorm2;
sol.J = cl_cost(opt,sol.uf(mm,:),sol.yf(pp,:))/Tv;
sol.Jd = (sol.J-sys.stoch_cost_Var)/opt.yrQTvnorm2;
if type == 6
    sol.J_yu = sol.J_yu / Tv;
    sol.J_reg = sol.J_reg / Tv;
end

end



% This stub function is invoked to launch all the DDPC approaches and save
% their step-by-step solution in the struct prd
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given data set contained in dpc
% - the initial condition for the state xini at each step
% - the type of predictive control approach
% - the prediction parameters and references contained in prd


function [prd] = ol(clx,dpc,xini,type,prd,beta,t)

prd.beta = 0;

switch type
    case 0
        prd = ol_GT(clx,dpc,prd,t);
    case 1
        prd = ol_MPC(clx,xini,prd,[],t);
    case 2
        prd = ol_2(clx,dpc,prd,beta,t);
    case 3
        prd = ol_3(clx,dpc,prd,beta,t);
    case 4
        prd = ol_23(clx,dpc,prd,beta,t);
    case 5
        prd = ol_DeePC(clx,dpc,prd,t);
    case 6
        prd = ol_FCE(clx,dpc,prd,t);
    case 7
        prd = ol_thm3(clx,dpc,prd,t);
end

end
