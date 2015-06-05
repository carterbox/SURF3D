function [maximums] = nonmaximumsupression(A, nbhood)
% NONMAXIMUMSUPRESSION takes a multidimensional array, A, and returns the
% locations and the values of the local maximums within distance defined by
% the vector NBHOOD.
%
% INPUTS
% A: The MxNxOxP dimensional array.
% neighborhood: a vector [n1,n2,n3,...] that defines the minimum spacing
% between local maxima. Analagous to the radius not the diameter.
%
% OUTPUTS
% maximums: an MxNDIMS+1 array with the coordinates of the maximums as row
% vectors and with their values in the last column.
%
% NOTES
% http://en.wikipedia.org/wiki/Flip_clock
%% -----------------------------------------------------------------------

% Pad the array so that its dimensions are a multiples of 2n+1.
% TODO: throw a warning instead of asserting and pad neighborhood with 1s.
assert(ndims(A) == ndims(nbhood),...
    'The dimensions of NBHOOD does not match A.');
if(iscolumn(nbhood)), nbhood = nbhood'; end;
a0  = size(A);
r = rem(a0,2.*nbhood+1);
A = padarray(A,r,'post'); % Post pad to avoid reindexing at the end.
a = size(A);
assert(sum(rem(a,2.*nbhood+1)) == 0);

% Generate the indicies of the points spaced at 2n+1 along the dimensions
% of the volume.
p = cell(ndims(A),1);
for i = 1:ndims(A);
    p{i} = 1+nbhood(i):2*nbhood(i)+1:a(i);
end

% Create a grid from these spaced indecies.
[X{1},X{2},X{3},X{4}] = ndgrid(p{1},p{2},p{3},{4});
numnbhoods = numel(X{1});
for j = 1:4, X{i} = reshape(X{i}, [numnbhoods,1]); end

% Break A into uniformly spaced neighborhoods 2n+1^dim and find the location of the
% partial maximum in each neighborhoods. 




% Then make sure none of these partial maximums are within the critical
% 2n+1 of each other.




end

%% Helper Functions ------------------------------------------------------

function [pmax] = partialmaximum(region)
% PARTIALMAXIMUM
%
%% -----------------------------------------------------------------------

end

function [X] = generategrid(p)
% GENERATEGRID does the same thing that matlab's NDGRID does, but the input
% and output are compressed into a single cell allowing for more flexible
% implementations.
%
% INPUTS
% p: a cell of vectors of indecies along the size of the grid.
%
% OUTPUT
% X: an array of each of the points in the grid as row vectors. 
%
% NOTES
% The inspiration for this flexible grid machine is the flip clock. I made
% this function because I didn't know ndgrid existed, and I'm too proud
% of my work to delete it. My implementation works by using the number of
% indecies along each dimension as a clock that counts down in order to
% iterate through all the possible combinations of indecies.
% http://en.wikipedia.org/wiki/Flip_clock
%% -----------------------------------------------------------------------

% Get the number of indecies for each dimension and allocate the output
% array.
flipclock = cellfun(@length, p); flipclock0 = flipclock;
X = zeros(prod(flipclock), length(flipclock));

for i = 1:length(X)    
    % Construct each point from the current time on the clock.
    for k = 1:length(flipclock)
        X(i,k) = p{k}(flipclock(k));
    end
    % Increment the time. We don't keep track of when flipclock(1) reaches
    % zero because we already know the total number of timesteps.
    for j = 0:length(flipclock)-2
        flipclock(end-j) = flipclock(end-j) - 1;
        if(flipclock(end-j) == 0)
            flipclock(end-j) = flipclock0(end-j);
        else
            break; % Break if we aren't changing anymore digits.
        end
    end
end

end