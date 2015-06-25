function [descriptors, isvalid] = extractFeatures3D(V,points,scale)
%EXTRACTFEATURES3D extracts SURF-162 descriptors for each of the points in the
%   volume V.
%
% INPUTS
% V: a 3D greyscale volume.
% points: a Nx3 matrix with each of the interest points as a row vector.
% [x,y,z,]
% scale:
%
% OUTPUTS
% descriptors: an NxM matrix where each row is a descriptor of length M for
% the nth point.
% isvalid: a Nx1 boolean matrix telling whether the nth point has a
% calculatable descriptor.
%
% NOTES
% This work is an extension of the work by Bay, Herbert, et al. "Speeded-up
% robust features (SURF)." Computer vision and image understanding 110.3 
% (2008): 346-359.
%
%% -----------------------------------------------------------------------

if(size(points,2) ~= 3), error('Points not 3 dimensional'), end;
if(ndims(V) ~= 3), error('Volume not 3 dimensional'), end;
scale = round(scale,2);
points = int16(points);

kNUMPOINTS = int32(size(points,1));
isvalid = true(kNUMPOINTS,1);

% Grab oriented subvolumes around each point.
[region,filters,ranges,hashtable,isvalid] = ...
                                    grabregions(V,points,scale,isvalid,20);

% Calculate the number of subregions and length of descriptor.
kSPLIT = int16([3,3,3]); % Change this line for SURF-384.
kNUMSUBREGIONS = prod(kSPLIT,'native');
kDESCRIPTORLENGTH = 6*kNUMSUBREGIONS;

% Generate the descriptors for each point with a valid region.
descriptors(kDESCRIPTORLENGTH,kNUMPOINTS) = single(0);
parfor j = 1:kNUMPOINTS
if(isvalid(j)) 
    % Get the hash value.
    h0 = hashtable(int16(scale(j)*100));
    
    % Make an integral image of the volume.
    J = integralimage3D(region{j});
    
    % Divide the integral image into equal subregions.
    [subregion, subfilter] = divideintegralimage(J,filters{h0},kSPLIT,kNUMSUBREGIONS);
    
    % Generate a sample spacing grid for Haar wavelets.
    [X,Y,Z] = ndgrid(ranges(h0,:),ranges(h0,:),ranges(h0,:));

    % Create the box filters for the wavelets.
    haarsize = 2*int16(scale(j));
    wavelets = makehaarwavets(haarsize);
    
    % TODO: Figure out a more elegant way to integrate the gaussian
    % weighting of the haar wavelet repsponses. The gaussian is centered at
    % the middle of the region.
    
    % Tabulate the reponse from each subregion and concatenate them
    descriptors(:,j) = computeSURFsums(subregion, subfilter, kNUMSUBREGIONS, kDESCRIPTORLENGTH, wavelets, X,Y,Z);
end
end

% Flip the output to the desired orientation.
descriptors = descriptors';

assert(length(descriptors) == length(points));
assert(length(isvalid) == length(points));
end
