function [color] = markkeypoints1( V0, points, outdir )
% MARKKEYPOINTS1 marks the POINTS with blue dots in the volume V0. Saves
% the new volume stack to OUTDIR/markedfeatures0001.png. #parallel
%
% INPUTS
% V0: an uint8 greyscale 3 dimensional array.
% points: an Nx3 array of key points [x,y,z].
%
% OUTPUTS
% color: a cell/stack of uint8 color images. 
% 
% NOTES
%
%% -----------------------------------------------------------------------
 
% Make an array the same size as V0 and put a pixel at each point.
marks = zeros(size(V0),'uint8');
for k = 1:size(points,1)
    marks(points(k,1),points(k,2),points(k,3)) = 255;
end

% Remove pixels in V0 at points so there are only pure colors.
mask = uint8(marks == 0);
V0 = mask.*V0;

% Create a stack to hold the color image slices.
color = cell(size(V0,3),1);
% Make a zero slice for the unused RGB channel.
z = zeros([size(V0(:,:,1),1), size(V0(:,:,1),2)],'uint8');
% Assmemble each of the slices. V0 is red, points are blue.
parfor i = 1:size(V0,3)
    color{i} = cat(3,V0(:,:,i),z,marks(:,:,i));
end

imstacksave(color,outdir,'markedfeatures');
end