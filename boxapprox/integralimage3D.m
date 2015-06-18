function [J] = integralimage3D(V)
% INTEGRALIMAGE3D computes the integral image of a 3D volume. Which is an
% array whose indecies contain the sum of all the intensities in the
% original image that are above and to the left. This one is inclusive and
% 1 indexed.
% 
% INPUT
% V: the volume.
%
% OUTPUT
% J (double): the integral image of V.
%
% NOTES
% http://en.wikipedia.org/wiki/Summed_area_table
% We need the extra bit depth for large volumes.
% Example: 128*1024^3 ~ 1.3x10^11 > 2^32
% Also, Matlab processes double faster than int64.
%% -----------------------------------------------------------------------

% Pad the integral image with zeros to make the loop simpler (There's no
% need for special cases at the edges). Make the default value -1 for
% troubleshooting purposes.
[x0,y0,z0] = size(V);
%J = padarray(-ones(x0,y0,z0, 'double'), [1,1,1], 'pre');
J(x0+1,y0+1,z0+1) = double(0);

for i = 2:x0+1
for j = 2:y0+1
for k = 2:z0+1
    % Single pass integral image computation extrapolated to 3D from
    % Wikipedia.
    J(i,j,k) = V(i-1,j-1,k-1) + J(i-1,j,k) + J(i,j-1,k) + J(i,j,k-1)...
    - J(i,j-1,k-1) - J(i-1,j,k-1) - J(i-1,j-1,k) + J(i-1,j-1,k-1);
end
end
end
end
