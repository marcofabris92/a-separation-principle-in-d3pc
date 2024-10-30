% This AUXILIARY stub function is employed to invoke the computation of a 
% noise-free ground-truth solution (not shown in the simulations; 
% it can be used to understand limit cases representing noise-free 
% scenarios) and it checks whether the result obtained is consistent 
% according to the Willems' fundamental lemma
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the parameters related to the given dataset contained in actual_dpc

% Invoked by:
% - MAIN.m
% Invokes:
% - cl.m, to launch the the closed-loop and perform the computations of
%   the ground-truth (selecting ol_GT.m)
% - lqdec_part.m, to compute the partition of the Hankel matrices


function [sol_GT] = get_GT(clx,actual_dpc)

% noise in the dynamics is forcedly turned down
clx.sys.flag_noise = 0;

% Hankel matrices
dpc.Up = actual_dpc.UpGT;
dpc.Uf = actual_dpc.UfGT;
dpc.Yp = actual_dpc.YpGT;
dpc.Yf = actual_dpc.YfGT;

% dpc_ini parameters
dpc.Tini = actual_dpc.Tini;
dpc.vec_uini = actual_dpc.vec_uini;
dpc.vec_yini = actual_dpc.vec_yini;
dpc.xini = actual_dpc.xini;

% LQ decomposition (for data... but data are noise-free in GT)
dpc.p1 = actual_dpc.p1;
dpc.p2 = actual_dpc.p2;
dpc.p3 = actual_dpc.p3;
dpc = lqdec_part(dpc,1);

% computation in advance of some repeatedly-used pieces of the solution
clx.L33_inv = eye(clx.opt.pT); 
clx.W12 = blkdiag(clx.opt.TR12,clx.opt.TQ12);
clx.W_star = clx.W12*[dpc.L22; dpc.L32];

% ground truth solution: (using noise free data)
sol_GT = cl(clx,dpc,0,0,[]);
check = norm(sol_GT.yf_hat-sol_GT.yf,"fro"); % this value must be ~ 0

fprintf('************************************\n')
fprintf('GROUND TRUTH NOISE-FREE ORACLE\n')
if 0 && check > 1e-10 && clx.opt.is_unconstrained
    error('Oracle is computed incorrectly.')
end
fprintf(['Noise-free oracle check: ' num2str(check) '\n'])
fprintf(['sol_GT.J = ' num2str(sol_GT.J) '\n'])
fprintf('************************************\n')
        
end