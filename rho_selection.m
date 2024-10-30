% This stub function finds the best values for the order rho of a model to 
% be ideintified; it computes
% - rho_GT, i.e. rho for a noise-free version of the dataset
% - rho, i.e. the list of model orders for each dataset
% - sys (it updates it)
% - the vector nsvars containing noise variances of the ARX estimation
% Given
% - the current version of the system parameters and references (sys)
% - the data sets stored in dat
% - the receding horizon T
% - the adopted criterion IC to validate ARX models
% - the flag rho_search_j: if = 1 then compute the list rho performing the 
%   ARX estimation for each dataset; if = 0 nothing occurs
% - the flag rho_search: if = 1 then choose the same order for 
%   each dataset and set it equal to rho_GT (the latter option can be
%   used to avoid to look for the best rho, since producing that list is
%   a quite slow operation that occurs once during the initialization);
%   if = 0 nothing occurs
% - the number of Monte Carlo runs NMC

% Invoked by:
% - MAIN.m, as the model order information is needed
% Invokes:
% - ARX_order.m, to perform the ARX estimation and find the best order
% according to the criterion in IC


function [rho_GT,rho,sys,nsvars] = ...
    rho_selection(sys,dat,T,IC,rho_search_j,rho_search,NMC)

rho_GT = ARX_order(sys,dat,T,IC,0);
rho = rho_GT*ones(1,NMC);
nsvars = sys.R*ones(1,NMC);
sys.rho_star = rho_GT;
if rho_search_j % selects a list of best orders rho
    parfor j = 1:NMC
        [rho(j),~,~,nsvars(j)] = ARX_order(sys,dat,T,IC,j);
    end
    % auxiliary quantity: average rho
    sys.rho_star = round(mean(rho));
end
if rho_search  % selects common rho_star if requested
    rho = sys.rho_star*ones(1,NMC);
end

end