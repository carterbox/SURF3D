function [J] = integralimage3D(V)
% INTEGRALIMAGE3D computes the integral image of a 3D volume. Which is an array
% whose indecies contain the sum of all the intensities in the original image
% that are above and to the left. This one is inclusive and 1 indexed.

V = double(V); % We need the extra bit depth for large volumes.
               % Example: 128*1024^3 ~ 1.3x10^11 > 2^32
               % Also, Matlab processes double faster than int64.

[u1,u2,u3] = size(V);
% Pad the integral image with zeros to make the loop simpler (There's no
% need for special cases at the edges). Make the default value -1 for
% troubleshooting purposes.
J = padarray(-ones(u1,u2,u3, 'double'), [1,1,1], 'pre');
for i = 2:u1+1
for j = 2:u2+1
for k = 2:u3+1
    
    % Single pass integral image computation extrapolated to 3D from
    % http://en.wikipedia.org/wiki/Summed_area_table
    
    J(i,j,k) = V(i-1,j-1,k-1) + J(i-1,j,k) + J(i,j-1,k) + J(i,j,k-1)...
    - J(i,j-1,k-1) - J(i-1,j,k-1) - J(i-1,j-1,k) + J(i-1,j-1,k-1);

end
end
end
end
