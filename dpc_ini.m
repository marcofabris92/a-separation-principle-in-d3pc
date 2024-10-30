% This function
% - provides the first initialization of the dpc struct, containing initial
%   trajectories for input, output and state
% - it resets parameters used by a Kalman filter to compute run ol_MPC.m 
% - it determines the dimensions of certain matrices (see below)
% Given
% - the underlying LTI system parameters contained in sys
% - the data sets stored in dat
% - the receding horizon T
% - the selected model order rho_star

% Invoked by:
% - MAIN.m, to initialize the struct dpc while computing the noise-less
%   ground-truth or the solution for all the MPC and DDPC approaches
% - outerMC.m, to reset dpc at each MC experiment if requested
% Invokes:
% - dynamics.m, to produce the initial trajectories


function [dpc,sys] = dpc_ini(sys,dat,T,rho_star)
    
% computation of initial trajectories for input, output and state
dpc.Ndata = dat.Ndata;
dpc.Tini = rho_star;              % prefixed-trajectory length
uini = zeros(sys.m,dpc.Tini);                           % uini
dpc.vec_uini = reshape(uini,[sys.m*dpc.Tini,1]);
[yini,xini] = dynamics(sys,sys.x0,uini,[],[]);          % yini
dpc.vec_yini = reshape(yini,[sys.m*dpc.Tini,1]);
dpc.xini = xini(:,end);                                 % xini

% reset KF initialization for the closed-loop (see cl.m)
sys.x_ini__ini_1 = sys.x0;
sys.P_ini__ini_1 = zeros(sys.n);
    
% exploiting initial information to get a better initialization for the KF
% that is used to compute the classic MPC solution (pretending to know the
% true model), see also ol_MPC.m
for t = 1:dpc.Tini
    L_t = sys.P_ini__ini_1*sys.C'*...
        ((sys.C*sys.P_ini__ini_1*sys.C'+sys.R)^-1);
    I_L_tC = sys.In-L_t*sys.C;

    x_filt_t = sys.x_ini__ini_1 +...
        L_t*(yini(:,t)-sys.C*sys.x_ini__ini_1-sys.D*uini(:,t));  
    P_filt_t = I_L_tC*sys.P_ini__ini_1*I_L_tC' + L_t*sys.R*L_t';

    sys.x_ini__ini_1 = sys.F*x_filt_t + sys.Z*yini(:,t) + ...
        sys.B_bar*uini(:,t);
    sys.P_ini__ini_1 = sys.F*P_filt_t*sys.F' + sys.Q_tilde;
end

% determination of the size of Hankel matrices (see Hankel.m)
dpc.TTini = dpc.Tini+T;              % # of Hankel rows
dpc.N = dat.Ndata-dpc.TTini+1;       % # of Hankel columns

% partitioning sizes for the LQ decomposition (see lqdec_part.m)
dpc.p1 = 1 : (sys.m+sys.p)*dpc.Tini;
dpc.p2 = dpc.p1(end)+1 : dpc.p1(end)+sys.m*T;
dpc.p3 = dpc.p2(end)+1 : dpc.p2(end)+sys.p*T;
    
end