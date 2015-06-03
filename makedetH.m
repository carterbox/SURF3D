function [detHgrid] = makedetH(J, filtersize)
% MAKEDETH takes an integral image and returns and array of the
% absolute values of the determinants of the Hessian at each calculatable
% point for the filtersize.
%
% INPUTS
% J: the integral image of some volume.
% filtersize: the size of the filter used to approximate H.
%
% OUTPUTS
% detHgrid: the array of |det(H)|.
%
% NOTES
% We take the absolute value of the det(H) because in three dimensions
% saddle points can also be detected. [Knopp et al. "Hough Transforms and 3D
% SURF for robust trhee dimenstional classification" pg. 4]
%% -----------------------------------------------------------------------

% First, approximate the Hessian matricies for all calculatable points in
% the volume.
gridofhessians = surfhessian3D(J, filtersize);

% Setup output cell. J is 1 larger than V;
[x,y,z] = size(J);
x0 = x-1; y0 = y-1; z0 = z-1;
detHgrid = zeros(x0,y0,z0,'double');

parfor i = 1:numel(gridofhessians)
    % Points near the edges were not calculated. Default value of detHgrid
    % is already zero.
    if(~isempty(gridofhessians{i}))
        % Form the full Hessian matrix and calculate |det(H)|.
        D = gridofhessians{i}
        H = [D(1),D(4),D(6);...
             D(4),D(2),D(5);...
             D(6),D(5),D(3)];
        detHgrid(i) = abs(det(H));
    end
end

end