% This function computes part of the approximation (51b) for the term 
% r_\gamma(u_f), see Theorem 3 of this paper. Its value is stored in the 
% variable omega. For more details see also Lemma 3, the Proof of
% Theorem 3 and the Proof of Theorem 4. In particular, equations (66) and 
% (67) are exploited for these computations.
% Given 
% - the parameters related to the given dataset contained in dpc
% - the receding horizon T
% - the (block diagonal) weight matrix Qo, weighting the output tracking
%   error ||yr-yf||^2

% Invoked by:
% - outerMC.m, for omega to be set
% Invokes:
% - J_(), to compute the matrices Jk defined in (66)
% - E_gi_gj(), to compute the expected value of g(i)*g(j) in (67)


function [omega] = reg_weight(dpc,T,Qo)

H = [dpc.Up; dpc.Yp; dpc.Uf];
Linv = pinv([dpc.L11 zeros(size(dpc.L11,1),size(dpc.L22,2));
             dpc.L21 dpc.L22]);
omega = zeros(size(Linv,2));
N = dpc.N;
pT = size(dpc.L33,2);
LQoL_N2 = dpc.L33'*Qo*dpc.L33/N^2;
parfor k_ = -T:T
    k = abs(k_);
    omega = omega + (N-k)*trace(LQoL_N2*J_(k,pT))*E_gi_gj(H,Linv,k);
end

end



% This function computes each matrix Jk defined in (66), where:
% Jk is the shift matrix such that [Jk]_i,j = 0 for j-i \neq k and
% [Jk]i,j = 1 for j-i == k
% Given
% - the index k
% - the size p of Jk

% Invoked by:
% - reg_weight.m, to perform the computation of omega
% Invokes: none


function J = J_(k,pT)
    
J = zeros(pT);
parfor i = 1:pT
    for j = 1:pT
        if j-i == k
            J(i,j) = 1;
        end
    end
end
    
end



% This function computes the expected value of g(i)*g(j) in (67)
% Given
% - the partial Hankel matrix H composed by past input/output trajectories
%   and the future input trajectories
% - the inverse matrix Linv with respect to the blocks L11, L21 and L22
% - the index k

% Invoked by:
% - reg_weight.m, to perform the computation of the regularizer omega
% Invokes: none


function result = E_gi_gj(H,Linv,k)

result = zeros(size(Linv,2));
for l = k+1:size(H,2)
    result = result + H(:,l)*H(:,l-k)';
end
result = Linv*result*Linv';

end