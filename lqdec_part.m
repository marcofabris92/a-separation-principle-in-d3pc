% This function computes the LQ decomposition and stores the result in dpc
% Given
% - the Hankel matrices contained in dpc
% - a flag called type, to distinguish between the ground_truth call and
%   a standard call to solve a gamma-DDPC optimization problem

% Invoked by:
% - get_GT.m, to compute a piece of the noise-free solution
% - MAIN.m, when the noise-free solution is required to be computed to 
%   provide a ground-truth (just for keep the code consistent)
% - outerMC.m, to compute pieces of solutions for gamma-DDPC approaches
% Invokes: 
% - lqdec(), see below


function [dpc] = lqdec_part(dpc,type)

switch type
    case 0 % noise-free ground truth
        [Q,L] = lqdec([dpc.UpGT; dpc.YpGT; dpc.UfGT; dpc.YfGT]);
        dpc.LGT11 = L(dpc.p1,dpc.p1);
        dpc.pinv_LGT11 = pinv(dpc.LGT11);
        dpc.LGT21 = L(dpc.p2,dpc.p1);
        dpc.LGT22 = L(dpc.p2,dpc.p2);
        dpc.LGT31 = L(dpc.p3,dpc.p1);
        dpc.LGT32 = L(dpc.p3,dpc.p2);
        dpc.LGT33 = L(dpc.p3,dpc.p3);
        dpc.Q_inv = pinv(Q);
%         dpc.QGT1 = Q(dpc.p1,:);
%         dpc.QGT2 = Q(dpc.p2,:);
%         dpc.QGT3 = Q(dpc.p3,:);
    case 1 % noisy gamma-DDPC
        [Q,L] = lqdec([dpc.Up; dpc.Yp; dpc.Uf; dpc.Yf]);
        dpc.L11 = L(dpc.p1,dpc.p1);
        dpc.pinv_L11 = pinv(dpc.L11);
        dpc.L21 = L(dpc.p2,dpc.p1);
        dpc.L22 = L(dpc.p2,dpc.p2);
        dpc.L31 = L(dpc.p3,dpc.p1);
        dpc.L32 = L(dpc.p3,dpc.p2);
        dpc.L33 = L(dpc.p3,dpc.p3);
        dpc.Q_inv = pinv(Q);
%         dpc.Q1 = Q(dpc.p1,:);
%         dpc.Q2 = Q(dpc.p2,:);
%         dpc.Q3 = Q(dpc.p3,:);
end

end



% This function returns the LQ decomposition of a given matrix (Q)
% Given
% - the matrix Q to be decomposed

% Invoked by:
% - lqdec_part.m, to be executed
% Invokes: none


function [Q,L] = lqdec(Q)
    [Q,L] = qr(Q',0);
    L = L';
    Q = Q';
end