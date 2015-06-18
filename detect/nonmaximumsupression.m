function [maximums] = nonmaximumsupression(A, nbhood)
% NONMAXIMUMSUPRESSION takes a 4-dimensional array, A, and returns the
% locations and the values of the local maximums within distance defined by
% the vector NBHOOD. #parallel
%
% INPUTS
% A: The MxNxOxP dimensional array.
% neighborhood: a col vector [n1;n2;n3;n4] that defines the minimum spacing
% between local maxima. Analagous to the radius not the diameter.
%
% OUTPUTS
% maximums: an Mx5 array with the coordinates of the maximums as row
% vectors and with their values in the last column.
%
% NOTES
% 
%% -----------------------------------------------------------------------

% Pad the array so that its dimensions are a multiples of 2n+1.
% TODO: throw a warning instead of asserting and pad neighborhood with 1s.
assert(ndims(A) == length(nbhood),...
    'The dimensions of NBHOOD does not match A.');
if(isrow(nbhood)), nbhood = nbhood'; end;
arraysize0  = size(A)';
remainders = rem(arraysize0,2.*nbhood+1);
A = padarray(A,double(2.*nbhood+1-remainders),'post'); % Post pad to avoid reindexing at the end.
arraysize = size(A)';
assert(sum(rem(arraysize,2.*nbhood+1)) == 0);

% Generate the indicies of the points spaced at 2n+1 along the dimensions
% of the volume.
p = cell(ndims(A),1);
for i = uint8(1:ndims(A));
    p{i} = 1+nbhood(i):2*nbhood(i)+1:arraysize(i);
end

% Create a grid from these spaced indecies.
[X{1},X{2},X{3},X{4}] = ndgrid(p{1},p{2},p{3},p{4}); clear('p');
num_nbhoods = numel(X{1});
for i = 1:4, X{i} = reshape(X{i}, [num_nbhoods,1]); end
X_ = cat(2,X{1},X{2},X{3},X{4}); clear('X');

% Break A into uniformly spaced neighborhoods 2n+1^dim and find the location of the
% partial maximum in each neighborhoods.
% TODO: Explicitly divide A so it is not a broadcast variable.
pmax_list = partialmaximum(X_,A,nbhood);

% Then make sure none of these partial maximums are within the critical
% 2n+1 of each other.
maximums = findneighbors(pmax_list, nbhood);

end

%% Helper Functions ------------------------------------------------------

function [pmax_list] = partialmaximum(X,A,nbhood)
% PARTIALMAXIMUM calculates the partial maximum for NBHOOD around X in A.
%
% INPUTS
% X: an array of the centers of all the neighborhoods to calculate partial
% maxima for.
% A: the array of data where we are searching for peaks.
% nbhood: a row vector denoting the size of the neighborhoods. Radii not
% diameters.
%
% OUTPUTS
% pmax_list: a Nx5 array of partial maximums where the first 4 columns are
% coordinates and the last column is the magnitude of a partial maximum.  
%% -----------------------------------------------------------------------

% Preallocate an output array.
pmax_list = zeros(size(X,1),size(X,2)+1);

% Carve out the region around the points in X and find the max in each.
% Record its location and magnitude.
numnbhoods = uint32(size(X,1));
for k = 1:numnbhoods
    range = cell(size(X,2),1);
    for i = 1:size(X,2)
       range{i} = (X(k,i) - nbhood(i)):(X(k,i) + nbhood(i));
    end
    thisregion = A(range{1},range{2},range{3},range{4});
    
    index = zeros(4,1);
    [magnitude, index(1)] = max(thisregion(:));
    [index(1),index(2),index(3),index(4)] = ind2sub(size(thisregion),index(1));
    index = X(k,:)' - nbhood - 1 + index;
    
    pmax_list(k,:) = [index', magnitude];
end

end

function [maximums] = findneighbors(pmax_list, nbhood)
% FINDNEIGHBORS locates the partial maximums that are in the neighborhood
% of others and removes the ones that are not local maximums.
%
% INPUTS
% pmax_list: a Nx5 array of partial maximums where the first 4 columns are
% coordinates and the last column is the magnitude of the max.  
%
% OUTPUTS
% maximums: a sorted list of the local maximums with a minimum spacing of
% 2n+1.
% 
% NOTES
% We actually average the components of nbhood to get a radius instead of
% using bounding boxes that are anisotropic.
%% -----------------------------------------------------------------------

% Search a KDtree for neighbors
tree = KDTreeSearcher(pmax_list(:,1:4));
radius = 2*mean(nbhood)+1;
matches = rangesearch(tree, pmax_list(:,1:4),radius);

% Remove points with matches that are larger.
num_pmax = uint32(size(pmax_list,1));
pmax_mags = pmax_list(:,5);
is_peak = true(num_pmax,1);
parfor i = 1:num_pmax
   % Compare each point with it's matches. No matches have length 0, and we
   % don't want peaks at zero.
   if(pmax_mags(i) == 0)
       is_peak(i) = false;
   else
       for j = 1:length(matches{i})
           % Mark the point if you find a match with a higher value.
           if(pmax_mags(i) < pmax_mags(matches{i}(j)))
               is_peak(i) = false;
               break;
           end
       end
   end
end

% Sort by magnitude and get rid of nonpeaks.
num_peaks = sum(is_peak);
pmax_list(:,5) = pmax_list(:,5).*is_peak;
pmax_list = sortrows(pmax_list, -5);
maximums = pmax_list(1:num_peaks,:);
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