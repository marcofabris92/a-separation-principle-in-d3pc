% This function determines the MPC solution for the model-based scheme by 
% means of an oracle (full knowledge of the underlying dynamics of the 
% system is here assumed) by relying on the Kalman filter (KF); 
% it computes
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the current state x_t of the system
% - the prediction parameters and references contained in prd
% - the possibility to force an input u, which can be set to [] 
%   if not needed

% Invoked by: 
% - cl.m, to run a comparison with the classic MPC, as if the input u
%   computed by other approaches could be used to control the closed-loop
%   in an MPC fashion (i.e. knowing the information of the model and 
%   exploiting Kalman filtering)
% - ol() in cl.m, to be started
% Invokes:
% - cvx_sol.m, if optimization constraints are enforced
% - prd_uy.m, to compute prd.u and prd.y
% - z_star.m, to compute gamma2


function [prd] = ol_MPC(clx,x_t,prd,u,t)

opt = clx.opt;
sys = clx.sys;
if ~isempty(u)
    prd.u = u;
end

%% computing the optimal uf, yf for the Kalman-based oracle,
% or assigning uf = u for a KF comparison

prd.y = sys.Gamma*prd.x_t__t_1;
if opt.is_unconstrained
    if isempty(u)
        % see data_generation.m for the definition of Upsilon matrices
        prd.u = sys.Upsilon_u*prd.ur + sys.Upsilon_y*(prd.yr-prd.y);
    end
    prd.y = prd.y + sys.Hd*prd.u;
else
    if isempty(u)
        prd = cvx_sol(clx,[],prd,1);
    else
        prd = cvx_sol(clx,[],prd,-1);
    end
end

%% getting the new KF-based prediction for the state

u_t = prd.u(1:sys.m);
y_t = sys.C*x_t+sys.D*u_t+sys.vy(:,t);
L_t = prd.P_t__t_1*sys.C'*((sys.C*prd.P_t__t_1*sys.C'+sys.R)^-1);
I_L_tC = sys.In-L_t*sys.C;

x_filt_t = prd.x_t__t_1 + L_t*(y_t-sys.C*prd.x_t__t_1-sys.D*u_t);
P_filt_t = I_L_tC*prd.P_t__t_1*I_L_tC' + L_t*sys.R*L_t';

prd.x_t__t_1 = sys.F*x_filt_t + sys.Z*y_t + sys.B_bar*u_t;
prd.P_t__t_1 = sys.F*P_filt_t*sys.F' + sys.Q_tilde;


end