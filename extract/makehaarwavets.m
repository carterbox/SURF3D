function boxpositions = makehaarwavets(filtersize)
%MAKEHAARWAVELETS generates the Haar wavelet filters for 3 directions at
%   the appropriate scale
%
% INPUTS
% scale (int16): the side length of the Haar filter. Should be even.
%
% OUTPUT
% boxpositions (int16): a cell of box filters. Each 7xB filter is in a cell
% and each box is represented as a column vector: [multiplier; corner; size].
%
% NOTES
% Because matlab is 1 indexed, the center of the filter is at
% filtersize/2 + 1.
%
%% -----------------------------------------------------------------------
assert(mod(filtersize,2) == 0);
l0 = filtersize/2;

boxpositions = cell(3,1);
parfor i = 0:2
    % Define the volumes for the wavelet
    box0 = circshift([-l0;-l0;-l0],i,1);
    box1 = circshift([  0;-l0;-l0],i,1);
    boxsize = circshift([l0;filtersize;filtersize],i,1);
    
    % Save it to the output
    boxpositions{i+1} = [-1,        1;          % multiplier
                          box0,     box1;       % corner
                          boxsize,  boxsize];   % size
end
end