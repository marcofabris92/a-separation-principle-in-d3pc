% This function returns the number of dimensions of a multi-array
% Given
% - a multiarray

% Invoked by:
% - MAIN.m, to retrieve grid dimensions
% - track_err_fig.m, to retrieve dimensions of a given signal
% Invokes: none


function n_dimensions = Ls(multiarray)

n_dimensions = 0;
if ~isempty(multiarray)
    n_dimensions = length(size(multiarray));
end

end