README

% This package reproduces the Numerical Examples conducted in Section 7
% of the paper entitled:

% "Harnessing Uncertainty for a Separation Principle in Direct
% Data-Driven Predictive Control"

% Authors:
% Prof. Alessandro Chiuso* Univ. of Padua      <alessandro.chiuso@unipd.it>
% Dr.   Marco Fabris**     Univ. of Padua         <marco.fabris.1@unipd.it>
% Prof. Valentina Breschi  Eindhoven Univ. of Technology <v.breschi@tue.nl>
% Prof. Simone Formentin   Polit. di Milano    <simone.formentin@polimi.it>
%
% *  A. Chiuso is the main algorithm developer of these Matlab scripts.
% ** M. Fabris is the main software developer of these MatLab scripts.

% Research Article submitted to the journal "Automatica", ID: 23-1641
% Publication history:
% - first submission: December 23rd, 2023
% - second submission: August 2nd, 2024
% - accepted on: October 3rd, 2024
% - last submission: October 31st, 2024

% Abstract:
% *************************************************************************
% Model Predictive Control (MPC) is a powerful method for complex system 
% regulation, but its reliance on an accurate model poses many limitations 
% in real-world applications. Data-driven predictive control (DDPC) aims at
% overcoming this limitation, by relying on historical data to provide 
% information on the plant to be controlled. 
% In this work, we present a unified stochastic framework for direct DDPC, 
% where control actions are obtained by optimizing the Final Control Error 
% (FCE), which is directly computed from available data only and 
% automatically weighs the impact of uncertainty on the control objective.
% Our framework allows us to establish a separation principle for 
% Predictive Control, elucidating the role that predictive models and their
% uncertainty play in DDPC. Moreover, it generalizes existing DDPC methods,
% like regularized Data-enabled Predictive Control (DeePC) and 
% $\gamma$-DDPC, providing a path toward noise-tolerant data-based control 
% with rigorous optimality guarantees. The theoretical investigation is 
% complemented by a series of numerical case studies, revealing that the 
% proposed method consistently outperforms or, at worst, matches existing 
% techniques without requiring tuning regularization parameters as other 
% methods do. 
% *************************************************************************

% This script and the remaining functions are documented at the beginning
% of each code fragment.




% --- Description of the most relevant files of this package ---

% (1) MAIN.m
% This script starts the numerical simulations presented in this paper,
% based on a Monte Carlo experiments to test the compared control methods
% over closed-loop dynamics.
% Some parameters in it must be selected by the users.
% Important: to reproduce the same results of showed in Figures 1 and 2
% one has to set sys.kind_of_setup = 0 first (this action reproduces
% Setups 1 and 2) and then set sys.kind_of_setup = 1 (this action 
% reproduces Setup 3). 

% (2) model_selection.m
% This subroutine is used to select models, tracking references and 
% optimization parameters.
% Important: the output reference in this file has to be manually
% constructed in order to reproduce our numerical results. So, make sure to
% change as desired any time MAIN.m is launched.

% (3) data_generation.m
% This function has two purposes:
% - compute and adjust some system parameters;
% - generate synthetic data.

% (4) rho_selection.m
% This function selects the model order rho at each Monte Carlo run,
% according to some criterion (e.g. the corrected AIC).

% (5) outerMC.m
% This subroutine constitutes the main body of each Monte Carlo run for
% these simulations.

% (6) innerMC_grids.m
% This subroutine constitutes the main body of each Monte Carlo run for
% the offline grid optimization of regularization parameters.

% (7) cl.m
% This function executes the closed-loop dynamics.

% (8) dynamics.m, sys_noise.m
% These functions exectute the open-loop dynamics or single updates by
% possibily adding process and/or output noise.

% (9) ol_XXXXX.m / cvx_sol.m
% The files ol_XXXXX.m encode the compared predictive control strategies.
% All these methods can be categorized are data-drive predictive control
% techniques, except for ol_MPC.m. The latter is a Kalman-based oracle to 
% show the best performance pretending to know the "true" model.
% Alternatively, one may run cvx_sol.m to solve constrained optimization
% problems through CVX. 
% Caution: cvx_sol.m is still under development; so, do not fully rely on 
% it. It may contain bugs and it could be very inefficient in certain
% situations.

% (10) dpc_ini.m
% This function
% - provides the first initialization of the initial trajectories for 
%   input, output and state; 
% - it resets parameters used by a Kalman filter to compute run ol_MPC.m; 
% - it determines the dimensions for partitioning Hankel matrices in 
%   gamma-DDPC approaches.

% (11) pbsidopt_ini.m
% Check out this subroutine if you wish to dive deep in the computations
% of the proposed FCE approach. Here the regularized in (37c) is
% constructed step by step.

% (12) reg_weight.m
% This subroutine has a simular task with respect to pbsidopt_ini.m at 
% computing the regularization term according to the approximation 
% discussed in Theorem 3.

% (13) tuning.m, SD_2.m, SD_3.m
% These functions are used in the online tuning for beta2 and beta3 of 
% gamma-DDPC approaches. 

% (14) hist_rho.m, my_boxplot.m, weird_trajectories.m, track_err_fig3.m
% These pieces of code are used to depict Figures 1 and 2.

% (15) table_exe_time.m
% This script provides the values in Table 1.