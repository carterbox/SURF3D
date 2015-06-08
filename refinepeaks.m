function [refined] = refinepeaks(peaks, A, d, octave)
% REFINEPEAKS takes the coarse location of PEAKS in A and uses quadratic
% fitting function in order to approximate the true location of the peak
% along dimension D.
%
% INPUT
% peaks:
% A:
% d: the scalar denoting the direction along which to fit the quadratic
% function
%
% OUTPUT
%
% NOTES
% 
%% -----------------------------------------------------------------------

assert(d <= ndims(A));
results = zeros(size(peaks,1),2);
% X is the range over which we will fit a quadratic.
X = octave'; query_grid = X(1):0.5:X(end);
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

refined = cat(2,peaks(:,1:3),results);

end