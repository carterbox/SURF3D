function [ output_args ] = assignorientation( V, point )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% INPUTS
% point Nx4 row vectors [x,y,z,scale]

%% -----------------------------------------------------------------------



% calculate Haar wavelet responses in the x y and z directions in a radius
% of 6s with a sample spacing of s

% set wavelet size to 4s

% weight responses according to gaussian with sigma = 2s.

% Assign direction according to max sum of responses within a pi/3 sector.

end

