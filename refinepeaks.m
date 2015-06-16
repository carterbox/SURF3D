function [refined] = refinepeaks(peaks, A, ~, octave)
% REFINEPEAKS takes the coarse location of PEAKS in A and uses quadratic
% fitting function in order to approximate the true location of the peak
% along dimension D.
%
% INPUT
% peaks: the Mx5 array of peaks from A. Where each peak is the following
% row vector [x,y,z,scale,magnitude].
% A: the 4 dimensional volume from which the peaks were extracted.
% d: the scalar denoting the direction along which to fit the quadratic
% function.
% octave: a row vector containing the filter sizes in the octave.
%
% OUTPUT
% refined: the Mx5 array of peaks from A. Where the dimension d has been
% refined using a quadratic fitting function.
% 
% NOTES
% 
%% -----------------------------------------------------------------------
%assert(d <= ndims(A));

% Make an array to hold the new scale and magnitude values.
results = zeros(size(peaks,1),2);
X = octave'; % the range over which we will fit a quadratic.
query_grid = X(1):0.5:X(end); % select desired spacing to interpolate
for i = 1:size(peaks,1)
    % Extract the data from A and put it in a vector
    coords = peaks(i,1:3);
    Y = squeeze(A(coords(1),coords(2),coords(3),:));
    % Interpolate the between filter sizes to find the estimated maximum
    % response.
    fitted_curve = interpn(X,Y,query_grid,'cubic');
    %if(rem(500,i) == 0), figure, plot(X,Y,'o',query_grid,fitted_curve,'-'); end
    [a_max,a_index] = max(fitted_curve);
    % Save the results to the output variable.
    results(i,:) = [query_grid(a_index),a_max];
end

% Replace the old scale and magnitude with the new ones.
refined = cat(2,peaks(:,1:3),results);

end