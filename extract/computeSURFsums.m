function v = computeSURFsums(subregion, subfilter, kNUMSUBREGIONS, kDESCRIPTORLENGTH, wavelets, X,Y,Z)
%COMPUTERSURFSUMS sums the wavelet responses at points [X,Y,Z] using
%   the integral image J.
%% -----------------------------------------------------------------------

v(kDESCRIPTORLENGTH,1) = single(0);
for i = 1:kNUMSUBREGIONS
    v(i*6-5:i*6) = computeSURFsums_helper(subregion{i},subfilter{i},wavelets,X,Y,Z);
end
    
end

function v = computeSURFsums_helper(J, filter, wavelets, X,Y,Z)
%
% INPUTS
% J (double): the integral image of the volume to be sampled.
% filter (single): weighting coeffients for each point in the volume.
% wavelets (int16): the box positions of the filters to be applied at each
% of the points. A cell of box filters. Each 7xB filter is in a cell and 
% each box is represented as a column vector: [multiplier; corner; size]. 
% Assumes three filters with 2 boxes each. XYZ: the coordinates of the
% points to be sampled where [X(n),Y(n),Z(n)] is the coordinate of the nth
% point to be sampled.
%
% OUTPUT
% v (single): the summed responses of the filters over the region.
% [dx,dy,dz,|dx|,|dy|,|dz|]
%
%% -----------------------------------------------------------------------
% Calculate the six terms of the vector by summing the response at each
% point.
v(6,1) = single(0); % [dx,dy,dz,|dx|,|dy|,|dz|]
numpoints = length(X);
for k = 1:numpoints
    % Calculate the reponse at this point.
    for j = 1:3
        thispoint = filter(X(k),Y(k),Z(k)) * single(wavelets{j}(1,1).*sumintegralimage3D(wavelets{j}(2:4,1) + [X(k);Y(k);Z(k)],wavelets{j}(5:7,1), J) + ...
                                                    wavelets{j}(1,2).*sumintegralimage3D(wavelets{j}(2:4,2) + [X(k);Y(k);Z(k)],wavelets{j}(5:7,2), J));
        v(j) = v(j) + thispoint;
        v(j+3) = v(j+3) + abs(thispoint);
    end
end
end