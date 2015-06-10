function [gridofhessians] = surfhessian3D(J, filtersize)
% SURFHESSIAN3D approximates the hessian matrix by using box filters in 3D.
% Returns a grid of hessians for each calculatable point in the volume.
%
% INPUTS
% J: the integral image of the volume.
% filtersize: the integer side length of the cube filter used to
% approximate H.
% fspacing: space between points where the filter is applied. i.e. the
% inverse sampling rate.
% 
% OUTPUTS
% gridofhessians: a cell of hessians at calculatable points in the volume.
% Hessian matrix elements are in this order: D11, D22, D33, D12, D23, D13.
%
% NOTES
% A calculateable point is one that does not overlap the edges of the
% volume.
% [H Bay, A Ess, T Tuytelaars, and L Van Gool. "Speeded Up Robust Features"
% 2008]
%% -----------------------------------------------------------------------

% The filter size must be multiple of 3 and greater than 8 in size.
assert( mod(filtersize, 3) == 0 && filtersize >= 9 );

% Setup output cell. J is 1 larger than V;
[x,y,z] = size(J);
x0 = x-1; y0 = y-1; z0 = z-1; 

% Generate the boxpositions for each of the two types of filters.
filter11 = makefilter11(filtersize);
filter12 = makefilter12(filtersize);

% Generate a list of all the places to apply the filter.
% Ignore points around the edges where the filter will give bad results.
fspacing = 1;
buffer = (filtersize - 1)/2;
[X,Y,Z] = ndgrid(1+buffer:fspacing:x0-buffer,...
                   1+buffer:fspacing:y0-buffer,...
                   1+buffer:fspacing:z0-buffer);
tic
tempgrid = cell(numel(X),1);
parfor k = 1:numel(X)
    %tic
    center = [X(k),Y(k),Z(k)];
    % Calculate each of the 6 terms of the Hessian.
    H = zeros(6,1);
    for j = 1:3
        % TODO: Move this space rotation outside the loop or make 3
        % separate filters to speed up the loop by an order of magnitude.
        rotatedJ = shiftdim(J,j-1); % Rotate the space instead of the filter.

        i = 1; % Calculate the main diagonal of the Hessian.
        while (i < length(filter11))
            H(j) = H(j) + filter11{i}.*sumintegralimage3D(...
                          filter11{i+1} + center, filter11{i+2}, rotatedJ);
            i = i + 3;
        end
        h = 1; % Calculate the off terms of the Hessian.
        while (h < length(filter12))
            H(j+3) = H(j+3) + filter12{h}.*sumintegralimage3D(...
                              filter12{h+1} + center, filter12{h+2}, rotatedJ);
            h = h + 3;
        end
    end
    %toc
    tempgrid{k} = H;
    %gridofhessians{center(1), center(2), center(3)} = H;
end

% TODO: Figure out if we can just reshape tempgrid instead of copying the
% elements manually. SOLVED: You can if you use ndgrid instead of
% meshgrid. (Meshgrid swaps the u1 and u2 coordinates) and if (X,Y,Z)
% samples every point in the volume. In our case, it doesn't because
% fspacing can be != 1 and we neglect the points around the edges.
gridofhessians = cell(x0,y0,z0);
for k = 1:numel(X)
    gridofhessians{X(k),Y(k),Z(k)} = tempgrid{k};
end
toc
end

%% Helper Functions

% TODO: Fix these second partial derivatives so that they taper like the
% 3D gaussians do. i.e. make the z direction non constant.

% TODO: Add a parameter to rotate the filters using SHIFTDIM.

function [boxpositions] = makefilter11(filter_size)
% Returns the corner locations for the regions of the 2nd partial
% derivative box filter D11.
% 
% INPUTS
% filter_size: the scalar length of the sizes of the cubic filter
%
% OUTPUT
% boxpositions: a row cell with information about each of the boxes in
% this format {multiplier, corner, size, ...}
% ------------------------------------------------------------------------

% Determine the characterisitc length.
l0 = filter_size/3;

% Calculate the corner coordinates for each of boxes.
middle = [-(l0 - 1)/2, -(l0 - 1), -(l0 - 1)];
top = middle + [l0, 0, 0];
bottom = middle - [l0, 0, 0];

% Calcculate the size of each box.
boxsize = [l0, 2*l0 - 1, 2*l0 - 1];

% Assign multipliers and put into a cell.
boxpositions = {1, top,    boxsize,...
               -2, middle, boxsize,...
                1, bottom, boxsize};
             % {multiplier, corner, size, ...}
end

function [boxpositions] = makefilter12(filter_size)
% Returns the corner locations for the regions of the 2nd partial
% derivative box filter D12.
% 
% INPUTS
% filter_size: the scalar length of the sizes of the cubic filter
%
% OUTPUT
% boxpositions: a row cell with information about each of the boxes in
% this format {multiplier, corner, size, ...}
% ------------------------------------------------------------------------

% Determine the characterisitc length.
l0 = filter_size/3;

% Calcculate the size of each box.
boxsize = [l0, l0, 2*l0+1];

% Calculate the corner coordinates for each of boxes.
box0 = [-l0, 1,   -l0];
box1 = [-l0, -l0, -l0];
box2 = [ 1,  1,   -l0];
box3 = [ 1,  -l0, -l0];

% Assign multipliers and put into a cell.
boxpositions = {1, box0, boxsize,...
               -1, box1, boxsize,...
               -1, box2, boxsize,...
                1, box3, boxsize};
             % {multiplier, corner, size, ...}
end