function [V3] = markkeypoints( V0, points )
%MARKKEYPOINTS marks the POINTS with circles in the volume V0.
%   Points are marked with circles that are the radius of the scale of the
%   detected key point.
%
% INPUTS
% V0: an uint8 greyscale 3 dimensional array.
% points: an Nx4 array of key points [x,y,z,radius].
%
% OUTPUTS
% V2: a cell/stack of uint8 color images. 
% 
% NOTES
% The points are marked as circles existing in the 3 planes parallel to the
% principal directions not as spheres; spheres would be too confusing.
%
% The main for loop can't be done in parallel, but the drawing of the
% circle could be moved out of the while loop into a parfor by pre breaking
% the list of points into groups by slice.
%% -----------------------------------------------------------------------
tic
points = int32(points(1:floor(end/10),1:4));

% We have to convert the images to color in order to draw the circles so
% to do that each location needs a pixel depth of three. To do that we make
% a [MxNxO] cell array of [1x1x3] cells.
V1 = cell(size(V0));
z = zeros([size(V0(:,:,1),1), size(V0(:,:,1),2)],'uint8');
V1(:,:,1) = num2cell(cat(3,V0(:,:,1),z,z),3);
%V1(:,:,1) = num2cell(repmat(V0(:,:,1),[1,1,3]),3);
for i = 2:size(V0,3)
    V1(:,:,i) = num2cell(cat(3,V0(:,:,i),z,z),3);
    %V1(:,:,i) = num2cell(repmat(V0(:,:,i),[1,1,3]),3);
end

% Make a shape inserting system object.
color = uint8([0 0 255]); % It's blue.
shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',color);

V3 = cell(3,1);
for direction = 1:3
    % Sort the points so we can draw circles on the same slice at the
    % same time.
    points = sortrows(points,direction);
    V1_ = shiftdim(V1,1);
    
    start_point = 1;
    end_point = 2; % [start:end) the range of points we are going to draw.
    i_slice = 1; % The slice that we are trying to draw on.
    while(start_point <= size(points,1) && i_slice <= size(V1_,direction))
        % Determine if there are any circles on this slice.
        if(points(start_point,3) == i_slice)
            % Find all the points on this slice.
            while(end_point < size(points,1) && points(end_point,3) <= i_slice)
                end_point = end_point + 1;
            end
            this_points = points(start_point:end_point-1,1:4); 
                        
            % Convert slice into an RGB image.
            slice = cell2mat(V1_(:,:,i_slice));
            
            % Define the circles to be drawn.
            r = this_points(:,4);
            this_points = circshift(this_points(:,1:3),-direction,2);
            circles = [this_points(:,1:2) r];
            
            % Draw circles
            slice = step(shapeInserter, slice, circles);
        
            % Convert slice back into cell and copy it back to V1.
            V1_(:,:,i_slice) = num2cell(slice,3);
            
            start_point = end_point;
        end
        i_slice = i_slice + 1; 
    end
    
    % Convert to a stack of color images.
    V2 = cell(size(V1_,3),1);
    for i = 1:size(V1_,3)
        V2{i} = cell2mat(V1_(:,:,i));
    end
    
    V3{direction} = V2;
end


toc

imstacksave(V3{3},'./circles/3','divided');
imstacksave(V3{2},'./circles/2','divided');
imstacksave(V3{1},'./circles/1','divided');
end

