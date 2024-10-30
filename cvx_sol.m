% This is an AUXILIARY function that should be used to handle 
% constrained optimization problems. 
% It leverages CVX (see https://cvxr.com/cvx/) and it computes
% - the next T control inputs starting from time instant t (prd.u), 
%   from which only the first element of this array is employed to 
%   push forward the closed-loop
% - the corresponding next T output predictions (prd.y)
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given data set contained in dpc
% - the prediction parameters and references contained in prd
% - type: an integer that determines which appraoch should be used (see
%   below, at the beginning of cl()

% Invoked by:
% - ol_2.m
% - ol_3.m
% - ol_23.m
% - ol_DeePC.m
% - ol_FCE.m
% - ol_GT.m
% - ol_MPC.m
% - ol_thm3.m
% Invokes: none


function [prd] = cvx_sol(clx,dpc,prd,type)

global cvx_sol_time cvx_sol_counter
cvx_sol_counter = cvx_sol_counter + 1;
cvx_sol_time_tic = tic;


opt = clx.opt;
zmT = zeros(opt.mT,1);
zpT = zeros(opt.pT,1);

% cost functional
Ucost = "+(prd.ur-uf)'*opt.TR*(prd.ur-uf)";
Ycost = "+(prd.yr-yf_hat)'*opt.TQ*(prd.yr-yf_hat)";
fnct = strcat(Ucost,Ycost);
Ucost_DeePC = "+(prd.ur-dpc.Uf*g)'*opt.TR*(prd.ur-dpc.Uf*g)";
Ycost_DeePC = "+(prd.yr-dpc.Yf*g)'*opt.TQ*(prd.yr-dpc.Yf*g)";
fnct_DeePC = strcat(Ucost_DeePC,Ycost_DeePC);

% the three regularization added to the cost functional
reg2 = "+prd.beta*gamma2'*gamma2";
reg3 = "+prd.beta*gamma3'*gamma3";
gamma2z = "gamma2 == zmT;";
gamma3z = "gamma3 == zpT;";
if type == 4
    reg2 = "+prd.beta(1)*gamma2'*gamma2";
    reg3 = "+prd.beta(2)*gamma3'*gamma3";
end
reg1Pi = "+prd.lmb1*g'*dpc.IPi2*g+prd.lmb2*sum(q)";
zg = zeros(dpc.N,1);
gis0 = "g == zg;";
if prd.lmb1 == +Inf
    reg1Pi = "+prd.lmb2*sum(q)";
end

% ddpc constraints
data_uf = "uf == dpc.L21*prd.gamma1+dpc.L22*gamma2;";
data_yf_hat = "yf_hat == dpc.L31*prd.gamma1+dpc.L32*gamma2";
data_yf_hat33 = "+dpc.L33*gamma3";
%data_trajectories = "[prd.zini; uf; yf_hat] == dpc.H*g";
data_trajectories_past = "prd.zini == dpc.Zp*g";
if prd.lmb1 == +Inf
    data_trajectories_past = ...
        [data_trajectories_past "zeros(dpc.N,1) == dpc.Ipi*g"];
end

% KF-based oracle constraint
kf_io = "yf_hat == prd.y+clx.sys.Hd*uf;";

% 1-norm constraints
one_norm_plus = "g <= q;";
one_norm_minus = "-g <= q;";


% FCE
pbs = clx.pbs;
YcostFCE = "+(pbs.yc-pbs.PsiU_hat*uf)*pbs.Q*(pbs.yc-pbs.PsiU_hat*uf)";
yrzini = [prd.yr; prd.zini];
reg_FCE = "+clx.beta_0*[yrzini' uf']*pbs.LambdaDopt*[yrzini; uf]";
fnct_FCE = strcat(Ucost,YcostFCE);

% thm3
r1 = length(prd.gamma1);
Wthm3 = prd.Reg_Thm3(r1+1:end,r1+1:end);
reg_thm3 = "+gamma2'*Wthm3*gamma2";

switch type
    case -1 % invoked by ol_MPC.m (with preset input, i.e. u is not empty)
        uf = prd.u;
        cvx_begin quiet
        variables yf_hat(opt.pT)
        minimize(eval(fnct))
        subject to
        eval(kf_io)
        eval(opt.set_cnstr)
        cvx_end
    case 0 % invoked by ol_GT.m
        cvx_begin quiet
        variables gamma2(opt.mT) uf(opt.mT) yf_hat(opt.pT)
        minimize(eval(fnct))
        subject to
        eval(data_uf)
        eval(data_yf_hat);
        eval(opt.set_cnstr)
        cvx_end
        prd.gamma2 = gamma2;
    case 1 % invoked by ol_MPC.m
        cvx_begin quiet
        variables uf(opt.mT) yf_hat(opt.pT)
        minimize(eval(fnct))
        subject to
        eval(kf_io);
        eval(opt.set_cnstr)
        cvx_end
    case 2 % invoked by ol_2.m
        prd.gamma2 = zeros(opt.mT,1);
        if prd.beta < +Inf
            cvx_begin quiet
            variables gamma2(opt.mT) uf(opt.mT) yf_hat(opt.pT)
            minimize(eval(strcat(fnct,reg2)))
            subject to
            eval(data_uf)
            eval(data_yf_hat);
            eval(opt.set_cnstr)
            cvx_end
            prd.gamma2 = gamma2;
        else % feasibility check
            cvx_begin quiet
            variables gamma2(opt.mT) uf(opt.mT) yf_hat(opt.pT)
            minimize(eval(strcat(fnct,reg2)))
            subject to
            eval(data_uf)
            eval(data_yf_hat);
            eval(gamma2z);
            eval(opt.set_cnstr)
            cvx_end
        end
    case 3 % invoked by ol_3.m
        prd.gamma3 = zeros(opt.pT,1);
        if prd.beta < +Inf
            cvx_begin quiet
            variables gamma2(opt.mT) gamma3(opt.pT) uf(opt.mT) yf_hat(opt.pT)
            minimize(eval(strcat(fnct,reg3)))
            subject to
            eval(data_uf)
            eval(strcat(data_yf_hat,data_yf_hat33));
            eval(opt.set_cnstr)
            cvx_end
            prd.gamma2 = gamma2;
            prd.gamma3 = gamma3;
        else % feasibility check
            cvx_begin quiet
            variables gamma2(opt.mT) gamma3(opt.pT) uf(opt.mT) yf_hat(opt.pT)
            minimize(eval(strcat(fnct,reg3)))
            subject to
            eval(data_uf)
            eval(strcat(data_yf_hat,data_yf_hat33));
            eval(gamma3z);
            eval(opt.set_cnstr)
            cvx_end
            prd.gamma2 = gamma2;
        end
    case 4 % invoked by ol_23.m
        beta2 = prd.beta(1);
        beta3 = prd.beta(2);
        prd.gamma2 = zeros(opt.mT,1);
        prd.gamma3 = zeros(opt.pT,1);
        if beta2 < +Inf && beta3 < +Inf
            cvx_begin quiet
            variables gamma2(opt.mT) gamma3(opt.pT) uf(opt.mT) yf_hat(opt.pT)
            minimize(eval(strcat(fnct,reg2,reg3)))
            subject to
            eval(data_uf)
            eval(strcat(data_yf_hat,data_yf_hat33));
            eval(opt.set_cnstr)
            cvx_end
            prd.gamma2 = gamma2;
            prd.gamma3 = gamma3;
        end
        if beta2 < +Inf && beta3 == +Inf
            cvx_begin quiet
            variables gamma2(opt.mT) gamma3(opt.pT) uf(opt.mT) yf_hat(opt.pT)
            minimize(eval(strcat(fnct,reg2,reg3)))
            subject to
            eval(data_uf)
            eval(strcat(data_yf_hat,data_yf_hat33));
            eval(gamma3z)
            eval(opt.set_cnstr)
            cvx_end
            prd.gamma2 = gamma2;
        end
        if beta2 == +Inf && beta3 < +Inf
            cvx_begin quiet
            variables gamma2(opt.mT) gamma3(opt.pT) uf(opt.mT) yf_hat(opt.pT)
            minimize(eval(strcat(fnct,reg2,reg3)))
            subject to
            eval(data_uf)
            eval(strcat(data_yf_hat,data_yf_hat33));
            eval(gamma2z);
            eval(opt.set_cnstr)
            cvx_end
            prd.gamma3 = gamma3;
        end
        if beta2 == +Inf && beta3 == +Inf
            cvx_begin quiet
            variables gamma2(opt.mT) gamma3(opt.pT) uf(opt.mT) yf_hat(opt.pT)
            minimize(eval(strcat(fnct,reg2,reg3)))
            subject to
            eval(data_uf)
            eval(strcat(data_yf_hat,data_yf_hat33));
            eval(gamma2z);
            eval(gamma3z);
            eval(opt.set_cnstr)
            cvx_end
        end
    case 5 % invoked by ol_DeePC.m
        uf = zeros(opt.mT,1);
        yf_hat = zeros(opt.pT,1);
        if prd.lmb1 < +Inf
            cvx_begin quiet
            variables g(dpc.N) q(dpc.N)
            minimize(eval(strcat(fnct_DeePC,reg1Pi)))
            subject to
            eval(data_trajectories_past);
            eval(one_norm_plus)
            eval(one_norm_minus)
            eval(opt.set_cnstr)
            cvx_end
            uf = dpc.Uf*g;
            yf_hat = dpc.Yf*g;
        else % feasibility check
            cvx_begin quiet
            variables g(dpc.N) q(dpc.N)
            minimize(eval(strcat(fnct_DeePC,reg1Pi)))
            subject to
            eval(data_trajectories_past);
            eval(gis0);
            eval(one_norm_plus)
            eval(one_norm_minus)
            eval(opt.set_cnstr)
            cvx_end
        end
    case 6 % invoked by ol_FCE.m
        cvx_begin quiet
        variables uf(opt.mT)
        minimize(eval(strcat(fnct_FCE,reg_FCE)))
        subject to
        eval(opt.set_cnstr)
        cvx_end
        yf_hat = pbs.I_PsiY_hat_inv*(pbs.yp+pbs.PsiU_hat*uf);
    case 7 % invoked by ol_thm3.m
        cvx_begin quiet
        variables gamma2(opt.mT) uf(opt.mT) yf_hat(opt.pT)
        minimize(eval(strcat(fnct,reg_thm3)))
        subject to
        eval(opt.set_cnstr)
        cvx_end
end

% eval(opt.set_cnstr)

% cvx_slvitr
% cvx_status

prd.u = uf;
prd.y = yf_hat;

cvx_sol_time = cvx_sol_time + toc(cvx_sol_time_tic);

end