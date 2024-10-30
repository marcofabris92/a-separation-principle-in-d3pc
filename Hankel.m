% This is a stub function to construct Hankel matrices within the context 
% of this research; it returns
% - a Hankel matrix H, given a series
% - its partition Hp, with respect to a past trajectory of a given series
% - its partition Hf, with respect to a future trajectory of a given series
% Given
% - the parameters related to the given data set contained in dpc
% - a numerical sequence (series)

% Invoked by:
% - MAIN.m, to compute the noise-free ground truth
% - outerMC.m, to compute the Hankel matrices for the j-th dataset
% Invokes:
% - Hankel_construction(), see below


function [Hp,Hf,H] = Hankel(dpc,series)

H = Hankel_construction(series,dpc.TTini,dpc.N)/sqrt(dpc.N);
d_s = size(series,1);
Hp = H(1:d_s*dpc.Tini,:);
Hf = H(1+d_s*dpc.Tini:end,:);
    
end




% This function builds up a Hankel matrix H
% Given 
% - a time series
% - the number of rows N_rows of H
% - the number of columns N_cols of H

% Invoked by:
% - Hankel.m, to perform the actual computations of H
% Invokes: none


function H = Hankel_construction(series,N_rows,N_cols) 

d_s = size(series,1);
H = zeros(N_rows*d_s,N_cols);
for i = 1:N_rows
    ii = (1:d_s)+(i-1)*d_s;
    for j = 1:N_cols
        H(ii,j) = series(:,i+j-1);
    end
end
    
end