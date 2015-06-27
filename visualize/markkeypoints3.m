function [] = markkeypoints3( V0, points, outdir )
%MARKKEYPOINTS3 cuts subvolumes out of V0 around the POINTS and saves them
% as a single image to a directory OUTDIR/p000001.png. #parallel
%
% INPUTS
% V0: an uint8 greyscale 3 dimensional array.
% points: an Nx3 array of key points [x,y,z,diameter].
% 
% NOTES
% Subvolumes are saved as a single image by concatinating the slices from
% left to right.
%% -----------------------------------------------------------------------

pad = ceil((max(points(:,4))-1)/2);
padded_V = padarray(V0,[pad,pad,pad],'both');
padded_points = uint32(points(:,1:3) + pad);
r = (points(:,4)-1)/2;

parfor i = 1:size(points,1)
    
    xrange = padded_points(i,1)-r(i):padded_points(i,1)+r(i);
    yrange = padded_points(i,2)-r(i):padded_points(i,2)+r(i);
    zrange = padded_points(i,3)-r(i):padded_points(i,3)+r(i);
    subvol = padded_V(xrange,yrange,zrange);
    
    slide = subvol(:,:,1);
    for j = 2:size(subvol,3)
        slide = cat(2,slide,subvol(:,:,j));
    end
    
    % Save the volume to a named directory
    folder = sprintf('/p%05i.png',i);    
    imwrite(slide,[outdir folder]);

end
end