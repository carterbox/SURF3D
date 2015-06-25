function [gridofhessians] = surfhessian3D(J, filtersize)
%SURFHESSIAN3D approximates the hessian matrix by using box filters in 3D.
%   Returns a grid of hessians for each calculatable point in the volume.
%   #parallel
%
% INPUTS
% J (double): the integral image of the volume.
% filtersize: the integer side length of the cube filter used to
% approximate H.
% fspacing: space between points where the filter is applied. i.e. the
% inverse sampling rate.
% 
% OUTPUTS
% gridofhessians (single): a cell of hessians at calculatable points in the
% volume. Hessian matrix elements are in this order: D11, D22, D33, D12,
% D23, D13.
%
% NOTES
% A calculateable point is one that does not overlap the edges of the
% volume. Expect 0.0013 seconds per Hessian on a 1.4 GHz processor.
% [H Bay, A Ess, T Tuytelaars, and L Van Gool. "Speeded Up Robust Features"
% 2008]
%% -----------------------------------------------------------------------

% The filter size must be multiple of 3 and greater than 8 in size.
assert( mod(filtersize, 3) == 0 && filtersize >= 9 );

% Setup output cell. J is 1 larger than V;
[x,y,z] = size(J);
x0 = x-1; y0 = y-1; z0 = z-1; 

% Generate the boxpositions for each of the two types of filters.
filters = makefilters(filtersize);

% Generate a list of all the places to apply the filter.
% Ignore points around the edges where the filter will give bad results.
fspacing = int32(2);
buffer = int32((filtersize - 1)/2);
[X,Y,Z] = ndgrid(1+buffer:fspacing:x0-buffer,...
                 1+buffer:fspacing:y0-buffer,...
                 1+buffer:fspacing:z0-buffer);
             
% store these values in advance to prevent recalculation at each iteration
numpoints = uint32(numel(X));
tempgrid = cell(numpoints,1);
numfilters = uint8(cellfun(@size,filters,{2;2;2;2;2;2}));  

parfor k = 1:numpoints
    % Calculate each of the 6 terms of the Hessian.
    H = zeros(6,1,'single');
    for j = uint8(1):6
        for i = 1:numfilters(j)
            H(j) = H(j) + single(filters{j}(1,i).*sumintegralimage3D(...
              filters{j}(2:4,i) + [X(k);Y(k);Z(k)], filters{j}(5:7,i), J));
        end
    end
    tempgrid{k} = H;
    %gridofhessians{center(1), center(2), center(3)} = H;
end

% TODO: Figure out if we can just reshape tempgrid instead of copying the
% elements manually. SOLVED: You can if you use ndgrid instead of
% meshgrid. (Meshgrid swaps the u1 and u2 coordinates) and if (X,Y,Z)
% samples every point in the volume. In our case, it doesn't because
% fspacing can be != 1 and we neglect the points around the edges.
gridofhessians = cell(x0,y0,z0);
for k = 1:numpoints
    gridofhessians{X(k),Y(k),Z(k)} = tempgrid{k};
end
end

%% Helper Functions

function [boxpositions] = makefilters(filter_size)
% Returns the corner locations for the regions of the 2nd partial
% derivative box filter D11, D22, D33, D12, D23, D13.
% 
% INPUTS
% filter_size: the scalar length of the sizes of the cubic filter
%
% OUTPUT
% boxpositions (int32): a cell of box filters. Each 7xB filter is in a cell
% and each box is represented as a column vector: [multiplier; corner; size].
% ------------------------------------------------------------------------

% Determine the characterisitc length.
l0 = int32(filter_size/3);

boxpositions = cell(6,1);
for i = int8(0:2)
    % 11 22 33 Filters ---------------------------------------------------
    % Calculate the size of each box.
    boxsize = circshift([l0; 2*l0 - 1; 2*l0 - 1],i,1);
    
    % Calculate the corner coordinates for each of boxes.
    middle = circshift([-(l0 - 1)/2; -(l0 - 1); -(l0 - 1)],i,1);
    top = middle + circshift([l0; 0; 0],i,1);
    bottom = middle - circshift([l0; 0; 0],i,1);

    % Assign multipliers and put into a cell.
    boxpositions{i+1} = [1,      -2,       1;...      % multiplier
                         top,     middle,  bottom;... % corner
                         boxsize, boxsize, boxsize];  % size
                     
    % 12 23 31 filters ---------------------------------------------------
    % Calculate the size of each box.
    boxsize = circshift([l0; l0; 2*l0+1],i,1);

    % Calculate the corner coordinates for each of boxes.
    box0 = circshift([-l0; 1;   -l0],i,1);
    box1 =           [-l0; -l0; -l0];
    box2 = circshift([ 1;  1;   -l0],i,1);
    box3 = circshift([ 1;  -l0; -l0],i,1);

    % Assign multipliers and put into a cell.
    boxpositions{i+4} = [1,      -1,      -1,       1;...     % multiplier
                         box0,    box1,    box2,    box3;...  % corner
                         boxsize, boxsize, boxsize, boxsize]; % size
end
end
