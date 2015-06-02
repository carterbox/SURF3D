function [T] = sumintegralimage3D(corner,size,J)
% SUMINTEGRALIMAGE uses integral image J to calculate the integral over an
% area in an image defined by the coordinates CORNER and vector SIZE.

% Check for 3D vectors.
assert(length(corner) == 3 && length(size) == 3);

% Make sure the size of the area is not larger than J.
%assert()

% Define the corners of the cube in the integral image.
G = corner;
A = G + size;
B = G + [0 1 1].*size;
C = G + [0 0 1].*size;
D = G + [1 0 1].*size;
E = G + [1 1 0].*size;
F = G + [0 1 0].*size;
H = G + [1 0 0].*size;

% Add and subtract the corners according to this idea on Wikipedia
% http://en.wikipedia.org/wiki/Summed_area_table but extrapolated to 3D
T = J(A) - J(E) - J(D) - J(B) + J(C) + J(H) + J(F) - 3*J(G);

end
