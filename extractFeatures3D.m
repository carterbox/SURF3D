function [descriptors, isvalid] = extractFeatures3D(V,points)
% EXTRACTFEATURES3D extracts SURF descriptors for each of the points in the
% volume V.
%
% INPUTS
% V: a 3D greyscale volume.
% points: a Nx3 matrix with each of the interest points as a row vector.
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

















end