function [T] = sumintegralimage3D(corner,boxsize,J)
% SUMINTEGRALIMAGE3D uses integral image J to calculate the integral over
% an area in an image defined by the coordinates CORNER and vector BOXSIZE.
%
% INPUTS
% corner: the corner of the box closest to [0,0,0].
% boxsize: a vector containing the side lengths of the box.
% J: the integral image of some volume.
%
% OUTPUTS
% T: The sum of intensities of the points in the area of the represented
% volume.
%
% NOTES:
% http://en.wikipedia.org/wiki/Summed_area_table
% ------------------------------------------------------------------------

% Check for 3D vectors.
assert(length(corner) == 3 && length(boxsize) == 3);
% Make sure the size of the area is not larger than J.
assert(sum(sum(size(J) - boxsize)) >= 0);
% Convert to row vectors.
if(iscolumn(corner)), corner = corner'; end;
if(iscolumn(boxsize)), boxsize = boxsize'; end;

% Define the corners of the cube in the integral image.
G = corner;
A = G + boxsize;
B = G + [0 1 1].*boxsize;
C = G + [0 0 1].*boxsize;
D = G + [1 0 1].*boxsize;
E = G + [1 1 0].*boxsize;
F = G + [0 1 0].*boxsize;
H = G + [1 0 0].*boxsize;

% Add and subtract the corners according to this idea on Wikipedia
% http://en.wikipedia.org/wiki/Summed_area_table but extrapolated to 3D
% T = A - E - D - B + C + H + F - G
T = J(A(1),A(2),A(3))...
    - J(E(1),E(2),E(3)) - J(D(1),D(2),D(3)) - J(B(1),B(2),B(3))...
    + J(C(1),C(2),C(3)) + J(H(1),H(2),H(3)) + J(F(1),F(2),F(3))...
    - J(G(1),G(2),G(3));

end
