function [J] = integralImage3D(V)
% INTEGRALIMAGE3D computes the integral image of a 3D volume. Which is an array
% whose indecies contain the sum of all the intensities in the original image
% that are above and to the left. This one is inclusive and 1 indexed.

[u1,u2,u3] = size(V);
% pad the integral image with zeros to make the loop less complicated
J = zeros(size(V) + 1);

for i = 2:u1
for j = 2:u2
for k = 2:u3
    
    % Single pass integral image computation extrapolated to 3D from
    % http://en.wikipedia.org/wiki/Summed_area_table
    
    J(i,j,k) = V(i,j,k) + J(i-1,j,k) + J(i,j-1,k) + J(i,j,k-1)...
    - J(i,j-1,k-1) - J(i-1,j,k-1) - J(i-1,j-1,k) + J(i-1,j-1,k-1);

end
end
end

end
