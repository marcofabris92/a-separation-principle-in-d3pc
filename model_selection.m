% This subroutine is called to initialize the most relevant fields of the 
% structs sys and opt. It sets and computes several useful parameters and
% references, such as the state-space, the dimensions of variables,
% optimization weights, noise covariance matrices and some input-output 
% relations.
% Given
% - the struct sys to be updated
% - the closed-loop parameters and settings contained in clx
% - the struct opt to be updated
% - an index called model to select the preferred models among the 4 listed 
%   below... feel free to add another model here if you wish to test it

% Invoked by:
% - MAIN.m, to initalize the model parameters and references
% Invokes: none


function [sys,opt] = model_selection(sys,clx,opt,model)

sys.K = [];         % if non-empty: innovation model is chosen
sys.Lambda = [];    % if non-empty: innovation model is chosen

switch model
    case 0  % Dorfler - Coulson - Markovsky [15] (white noise output error)
        bf = [0        0       0  0.28261  0.50666];
        af = [1 -1.41833 1.58939 -1.31608  0.88642];
        [sys.A,sys.B,sys.C,sys.D] = tf2ss(bf,af);
        sys.Bn = zeros(4,1);
        sys.Dn = sqrt(sys.r);
    case 1  % Dorfler - Coulson - Markovsky [15] (colored noise inn. form)
%   ** this is the system chosen for the simulations shown in the paper **
        bf = [0        0       0  0.28261  0.50666];
        af = [1 -1.41833 1.58939 -1.31608  0.88642];
        [sys.A,sys.B,sys.C,sys.D] = tf2ss(bf,af);
%        figure
%        bode(tf(bf,af,1))
        sys.K = 1.9*[0.0939 -0.3433 0.1063 1.2058]';   
        sys.Lambda = 7.2880*1e-4;
    case 2  % Breschi - Chiuso - Formentin [8]
        sys.A = [0.7326 -0.0861;                  
                 0.1722  0.9909];
        sys.B = [0.0609; 
                 0.0064];
        sys.C = [0 1.4142];
        sys.D = 0;
        sys.Bn = zeros(2,1);
        sys.Dn = sqrt(sys.r);
    case 3 % Berberich [5]
        sys.A = [0.921 0     0.041     0;
                 0     0.918 0     0.033;
                 0     0     0.924     0; 
                 0     0     0     0.937];
        sys.B = [0.017 0.001;  
                 0.001 0.023;
                 0     0.061;
                 0.072     0];
        sys.C = [1 0 0 0;  
                 0 1 0 0];
        sys.D = zeros(2);
        sys.Bn = zeros(4,2);
        sys.Dn = sqrt(sys.r)*eye(2);
end

% dimensions and secondary paramters
sys.n = size(sys.A,1);      % number of states
sys.m = size(sys.B,2);      % number of inputs
sys.p = size(sys.C,1);      % number of outputs
sys.Ts = 1;                 % sampling time
sys.x0 = zeros(sys.n,1);    % initial condition (rewritten in data_gen.)
sys.In = eye(sys.n);
sys.na = ones(sys.p,sys.p); % A structure for arx model
sys.nb = ones(sys.p,sys.m); % B structure for arx model
sys.nk = ones(sys.p,sys.m); % delay structure for arx model

% optimization dimensions
opt.mT = sys.m*opt.T;
opt.pT = sys.p*opt.T;

% optimization weights
opt.R = opt.Wu*eye(sys.m);
opt.Q = opt.Wy*eye(sys.p);
opt.TR = kron(opt.I_T,opt.R);
opt.TQ = kron(opt.I_T,opt.Q);
opt.TR12 = chol(opt.TR);
opt.TQ12 = chol(opt.TQ);
opt.W12 = blkdiag(opt.TR12,opt.TQ12);

% references
sys.TvT = clx.Tv+opt.T-1;                 	   % reference dimension
opt.ur = zeros(sys.m,sys.TvT);                 % input reference
%opt.yr = sin((5*pi/(sys.TvT-1))*(1:sys.TvT)); % an old output reference
                                               % used in [7] and [9]

% square wave output reference used in the experiments of this paper
% for setup 1 and 2:
opt.yr = square((2*pi/40)*(1:sys.TvT));
opt.yr = repelem(opt.yr,sys.p,1);

% "irregular" square wave (output reference)
% CAVEAT: this input does not adapt to any Tv, it is assigned only for 
% Tv = 500
if sys.kind_of_setup % <-- 1: replicates setup 3
    opt.yr = zeros(1,sys.TvT);
    for t = 1:10
        opt.yr(t) = -1;
    end
    for t = 11:20
        opt.yr(t) = 10;
    end
    for t = 21:30
        opt.yr(t) = -1;
    end
    for t = 31:40
        opt.yr(t) = 1;
    end
    for t = 41:180
        opt.yr(t) = -10;
    end
    for t = 181:320
        opt.yr(t) = 10;
    end
    for t = 321:460
        opt.yr(t) = -1;
    end
    for t = 461:470
        opt.yr(t) = 1;
    end
    for t = 471:480
        opt.yr(t) = -10;
    end
    for t = 481:490
        opt.yr(t) = 1;
    end
    for t = 491:500
        opt.yr(t) = -1;
    end
end

% remove this if you wish to assign the "irregular" square wave
if sys.kind_of_setup
    opt.yr = sin((5*pi/(sys.TvT-1))*(1:sys.TvT));
end


opt.yrTvnorm2 = norm(opt.yr(:,1:clx.Tv),"fro")^2; % norm^2 of y_r (over Tv)
opt.yrQTvnorm2 = 0;
for t = 1:clx.Tv
    opt.yrQTvnorm2 = opt.yrQTvnorm2 + opt.yr(:,t)'*opt.Q*opt.yr(:,t);
end

% stacked references
opt.vec_ur = reshape(opt.ur,[sys.m*sys.TvT,1]);
opt.vec_yr = reshape(opt.yr,[sys.p*sys.TvT,1]);

% noise covariance initialization     
if isempty(sys.K)
%     uncomment and replace if Bn and Dn are not defined
%     sys.R = sys.r*eye(sys.p);   
%     sys.Q = sys.q*eye(sys.n);   
%     sys.S = zeros(sys.n,sys.p); 

    sys.R = sys.Dn*sys.Dn';  % output noise var., must be symm and pd
    sys.Q = sys.Bn*sys.Bn';  % state noise var., must be symm and psd
    sys.S = sys.Bn*sys.Dn';  % cross covariance
else
    sys.R = sys.Lambda;
    sys.S = sys.K*sys.R;
    sys.Q = sys.S*sys.K';
end


% input-output relations (Gamma and Hd), lag of the system (see [9])
sys.Gamma = zeros(opt.pT,sys.n);
sys.Hd = zeros(opt.pT,opt.mT);
ii = 1:sys.p;
jj = 1:sys.m;
sys.Gamma(ii,:) = sys.C;
sys.Hd(ii,jj) = sys.D;
sys.lag = 0;
for t = 1:opt.T-1
    new_entry = sys.Gamma(ii,:)*sys.A;
    ii = ii + sys.p;
    jj = jj + sys.m; 
    sys.Gamma(ii,:) = new_entry;
    sys.Hd(ii,jj) = sys.D;
    if sys.lag == 0 && rank(sys.Gamma(1:ii(1)-1,:)) == sys.n
        sys.lag = t;
    end
end
for low_diag = 0:opt.T-2
    ii = (1:sys.p)+low_diag*sys.p;
    new_entry = sys.Gamma(ii,:)*sys.B;
    for col = low_diag:opt.T-2
        ii = (1:sys.p)+(col+1)*sys.p;
        jj = (1:sys.m)+(col-low_diag)*sys.m;
        sys.Hd(ii,jj) = new_entry;
    end
end

end
