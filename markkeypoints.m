function [V1] = markkeypoints( V0, points )
%MARKKEYPOINTS marks the POINTS with circles in the volume V0.
%   Points are marked with circles that are the radius of the scale of the
%   detected key point. They are marked as circle existing in the 3 planes
%   parallel to the principal directions not as spheres. That would be too
%   confusing.
% INPUTS
%
% OUTPUTS
% 
%% -----------------------------------------------------------------------

% We have to convert the images to color in order to draw the circles so
% to do that each location needs a pixel depth of three. To do that we make
% a [MxNxO] cell array of [1x1x3] cells.
V1 = cell(size(V0));
V1(:,:,1) = num2cell(repmat(V0(:,:,1),[1,1,3]),3);
for i = 2:size(V0,3)
    V1(:,:,i) = num2cell(repmat(V0(:,:,i),[1,1,3]),3);
end

for direction = 1:3
    % Sort the points so we can draw circles on the same slice at the
    % same time.
    points = sortrows(points,direction);
    V1 = shiftdim(V1,1);
    
    i_points = 1; i_slice = 1;
    while(i_points <= size(points,1) && i_slice <= size(V1,direction))
        
        % Make a list of all the points on this slice.
        
        % Determine if there are any circles on this slice.
            % Convert this slice into an RGB image.
        
        [x,y,z,r] = points(i,:);
        
        shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',yellow);
        circles = [x y r];
        J = step(shapeInserter, slice, circles);
        
        V0(x,y,:) = V0(x,y,:);
    end
end

end

