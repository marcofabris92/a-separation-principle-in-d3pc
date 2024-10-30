% This function estimates the quantities Hs and (Hs*Lambda^(1/2))^-1, where
% Hs is defined in (7d) of [9], and Lambda is the covariance matrix of the
% innovation e(t), as defined in (3) of [9]. This quantity is then used
% in SD_2.m and SD_3.m for the tuning of beta2 and beta3, improving the
% methods proposed in [9].
% Given
% - the underlying LTI system parameters contained in sys
% - the data sets stored in dat
% - the optimization parameters and references contained in opt
% - the already selected model order
% - the number N of data to be selected from the j-th data set
% - the index j corresponding to the j-th Monte Carlo run

% Invoked by:
% - outerMC.m, to get ready before several calls to the different DDPC
% approaches
% Invokes: none


function [sys] = HsLam_estimation(sys,dat,opt,rho,N,j)

est_model = arx(iddata(dat.y_noisy(:,1:N,j)',...
    dat.u_noisy(:,1:N,j)'),[rho*sys.na rho*sys.nb sys.nk]);
sys.HsLam_hat_inv = zeros(opt.pT);
for t = 0:min(opt.T-1,rho)
    for r = 1:sys.p
        for c = 1:sys.p
            if sys.p > 1
                a_rc_t = est_model.a{r,c}(t+1);
            else
                a_rc_t = est_model.a(t+1);
            end
            for i = t+1:opt.T
                ii = (i-1)*sys.p+r;
                jj = (i-1-t)*sys.p+c;
                sys.HsLam_hat_inv(ii,jj) = a_rc_t;
            end
        end
    end
end
sys.Hs_hat = pinv(sys.HsLam_hat_inv);
sys.HsLam_hat_inv = pinv(kron(opt.I_T,sqrtm(est_model.noisevariance))) ...
    * sys.HsLam_hat_inv;

end