% This subroutine initializes the struct called pbs, which contains 
% all the parameters that are involved in the solution of the FCE DDPC
% approach (see also ol_FCE.m). Such computations rely on the 
% Predictor-Based Subspace identification (PBSid) method firstly devised 
% and analyzed in [11]. 
% It computes all the quantities needed to construct the regularizer r(uf) 
% provided in (37c), Theorem 2. To this purpose, we would like to provide
% more details about r(uf). In particular, we want to highlight explicitly
% the dependencies of r(uf,zini,yr) also on the initial input-output
% trajectories zini and the output reference yr.
% ************************************************************************
%
%                                   [Cyy Cyz | Cyu]   [ yr ]
% r(uf,zini,yr) = [yr' zini' uf'] * [Czy Czz | Czu] * [zini] 
%               |  \_______/        [Cuy Cuz | Cuu]   [ uf ]
%               |     =: v'
%               |
%               = [v' uf'] * [Cvv | Cvu] * [ v  ]
%                            [Cuv | Cuu]   [ uf ]
%                            \_________/
%                                =: C
% 
% Let's understand how to rewrite r(uf,zini,yr) into a quadratic form where
% it is emphasized the desired control uf* towards which the input uf is
% "pushed" by means of this kind of cost minimization. 
% 
% Start by putting C into a block diagonal form. Find a suitable
% similarity matrix:
% 
% [I T'] * [Cvv Cvu] * [I 0] = [Gvv   0]
% [0 I ]   [Cuv Cuu]   [T I]   [0   Cuu]
% 
% Cvu + T'*Cuu = 0  =>  T = -(Cuu^-1) * Cvz
%                                       \_/
%                                        = Czv'
% 
% [Cvv-Cvu*(Cuu^-1)*Cuv   0] * [I 0] = [Gvv   0]
% [Cuv                  Cuu]   [T I]   [0   Cuu]
% 
% Then write the cost as follows:
% 
% J(v,uf) = [v' uf'] * [Cvv Cvu] * [ v  ]
%         |            [Cuv Cuu]   [ uf ]
%         |
%         = [v' uf'] * [I -T'] * [I T'] * C * [I 0] * [I  0] * [ v  ]
%         |            [0  I ]   [0 I ]       [T I]   [-T I]   [ uf ]
%         |                                           \____________/
%         |                                               = [  v  ]
%         |                                                 [uf-Tv]
%         |
%         = [v' (uf-uf*)'] * [Gvv   0] * [  v   ]
%                            [0   Cuu]   [uf-uf*]
% 
% Thus, it must be
%                     uf* := Tv
%                          |
%                          = -(Cuu^-1)*Cuv*v
%                          |
%                          = -[(Cuu^-1)*Cuv*zini+(Cuu^-1)*Cuv*yr]
% 
% Finally, the cost to be minimized can be rewritten into the quadratic
% form yielded by
%                   J(uf,zini,yr) = v'*Gvv*v + (uf-uf*)'*Cuu*(uf-uf*)
%                                              \____________________/
%                                              only this part matters
%                                              when minimizing over uf
%
% ************************************************************************
% Given
% - the closed-loop variables, parameters and settings contained in clx
% - the data sets stored in dat
% - the selected model order rho
% - the index j_ addressing the j_-the Monte Carlo run

% Invoked by:
% - outerMC.m, to start this initialization
% Invokes: none

function pbs = pbsidopt_ini(clx,dat,rho,j_)

Ndata = dat.Ndata;
m = clx.sys.m;
p = clx.sys.p;
pm = p+m;
T = clx.opt.T;

%% upper bound equality constraint based on chi^2 stats
pbs.k1_a = p*T+1.96*sqrt(2*p*T);

%% estimation of theta hat and its variance
uD = dat.u_noisy(:,:,j_);
yD = dat.y_noisy(:,:,j_);

% racall that z := [y' u']'
ZD = zeros(Ndata-rho,rho*pm);
YD = zeros((Ndata-rho)*p,1);
for i = 1:Ndata-rho
    ii = (i-1)*p+1 : i*p;
    YD(ii) = yD(:,rho+i); 
    for j = 0:rho-1
        jj = pm*j+1 : pm*(j+1);
        ZD(i,jj) = kron(eye(p),[yD(:,rho+i-1-j)' uD(:,rho+i-1-j)']);
    end
end

pbs.sigma2_hat = clx.nsvars(j_);
theta_hat = pinv(ZD)*YD;
[~,SZD,VZD] = svd(ZD,"econ");
Sigma_hat_theta = pbs.sigma2_hat*VZD*pinv(SZD^2)*VZD';
Sigma_hat_theta = (Sigma_hat_theta+Sigma_hat_theta')/2;


%% permutation matrices

% P0
P0 = zeros(rho*p*pm);
counter = 1;
for k = 0:1
    Jk = p;
    if k
        Jk = m;
    end
    for j = 1:Jk
        for l = 1:rho
            for i = 1:p
                rn = rho*pm*(i-1)+pm*(l-1)+p*k+j;
                P0(counter,rn) = 1;
                counter = counter + 1;
            end
        end
    end
end
P0 = sparse(P0);

% N0
N0 = 1;
if rho > T
    N0 = sparse([eye(T) zeros(T,rho-T)]);
elseif rho < T
    N0 = sparse([eye(rho); zeros(T-rho,rho)]);
end

% Pz
Pz = zeros(rho*T*pm,rho*pm);
Nplus = zeros(rho);
for i = 1:rho-1
    Nplus(i,i+1) = 1;
end
for i = 1:rho
    ii = T*pm*(i-1)+1 : T*pm*i;
    Pz(ii,:) = kron(eye(pm),N0*(Nplus^(rho-i)));
end
Pz = sparse(kron(Pz,eye(p))*P0);


% Pu and Py
Py = zeros(p*T^2,rho*pm);
Pu = zeros(m*T^2,rho*pm);
for i = 1:T
    Nminus_i = zeros(T,rho);
    j = 1;
    while i+j <= T && j <= rho
        Nminus_i(i+j,j) = 1;
        j = j + 1;
    end
    ii = T*p*(i-1)+1 : T*p*i;
    Py(ii,1:rho*p) = kron(eye(p),Nminus_i);
    ii = T*m*(i-1)+1 : T*m*i;
    Pu(ii,rho*p+1:end) = kron(eye(m),Nminus_i);
end
Py = sparse(kron(Py,eye(p))*P0);
Pu = sparse(kron(Pu,eye(p))*P0);


%% computation of Xi, PsiU, I-PsiY and redefinition of Q
pbs.Xi_hat = reshape(Pz*theta_hat,[T*p,rho*pm]);
pbs.PsiU_hat = reshape(Pu*theta_hat,[T*p,T*m]);
pbs.PsiY_hat = reshape(Py*theta_hat,[T*p,T*p]);
pbs.I_PsiY_hat = eye(T*p)-reshape(Py*theta_hat,[T*p,T*p]);
pbs.I_PsiY_hat_inv = pbs.I_PsiY_hat^-1;

pbs.Q = pbs.I_PsiY_hat_inv'*clx.opt.TQ*pbs.I_PsiY_hat_inv;
J = (clx.opt.TQ12*pbs.I_PsiY_hat_inv)'; % J * J' = Q


%% computation of Lambda_D^opt
% Some of the following computations are skipped based on the
% technicalities discussed above in the documentation
% Bzz = Pz*Sigma_hat_theta*Pz';
% Bzu = Pz*Sigma_hat_theta*Pu';
% Bzy = Pz*Sigma_hat_theta*Py';
Buz = Pu*Sigma_hat_theta*Pz';
Buu = Pu*Sigma_hat_theta*Pu';
Buy = Pu*Sigma_hat_theta*Py';
% Byy = Py*Sigma_hat_theta*Py';

rz = rho*pm;
ru = T*m;
ry = T*p;

% Czz = zeros(rz);
% pbs.Czu = zeros(rz,ru);
% Czy = zeros(rz,ry);
% for i = 1:rz
%    ii = T*p*(i-1)+1 : T*p*i;
%     for j = 1:rz
%         jj = T*p*(j-1)+1 : T*p*j;
%         Czz(i,j) = trace(J'*Bzz(ii,jj)*J);
%     end
%     for j = 1:ru
%         jj = T*p*(j-1)+1 : T*p*j;
%         pbs.Czu(i,j) = trace(J'*Bzu(ii,jj)*J);
%     end
%     for j = 1:ry
%         jj = T*p*(j-1)+1 : T*p*j;
%         Czy(i,j) = trace(J'*Bzy(ii,jj)*J);
%     end
% end


pbs.Cuu = zeros(ru);
pbs.Cuy = zeros(ru,ry);
pbs.Cuz = zeros(ru,rz);
for i = 1:ru
    ii = T*p*(i-1)+1 : T*p*i;
    for j = 1:rz
        jj = T*p*(j-1)+1 : T*p*j;
        pbs.Cuz(i,j) = trace(J'*Buz(ii,jj)*J);
    end
    for j = 1:ru
        jj = T*p*(j-1)+1 : T*p*j;
        pbs.Cuu(i,j) = trace(J'*Buu(ii,jj)*J);
    end
    for j = 1:ry
        jj = T*p*(j-1)+1 : T*p*j;
        pbs.Cuy(i,j) = trace(J'*Buy(ii,jj)*J);
    end
end
pbs.Cuu = (pbs.Cuu+pbs.Cuu')/2;

% Cyy = zeros(ry);
% for i = 1:ry
%     ii = T*p*(i-1)+1 : T*p*i;
%     for j = 1:ry
%         Cyy(i,j) = trace(J'*Byy(ii,jj)*J);
%     end
% end

pbs.LambdaDopt = [zeros(rz)    pbs.Cuz' zeros(rz,ry);
                  pbs.Cuz      pbs.Cuu  pbs.Cuy;
                  zeros(ry,rz) pbs.Cuy'  zeros(ry)];


end