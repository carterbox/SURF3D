
function [region,filters,ranges,hashtable,isvalid] = grabregions(V,points,scale,isvalid,windowsize)
%GRABREGIONS segments oriented regions from the Volume around each point.
% It returns the regions, a table of filters and sampling ranges for each
% region size, and an array telling whether each descriptor is valid.
%
% INPUTS
% V (3D Array): The volume from which the regions come.
% points (int16): Mx3 matrix of coordinates of points in the Volume.
% scale (): the scale of each point
% isvalid (logical): an Mx1 vector telling whether each point has a valid
% descriptor
% windowsize: size of each region to be grabbed is WINDOWSIZE*SCALE + 1.
% WINDOWSIZE should be even.
%
% OUTPUTS
% region (cell): a cell of the segmented regions.
% filters (cell): a cell of gaussian filters; one for each unique scale.
% ranges (uint8): an Nx5 table of subregion sample spacings for the Haar
% wavelet responses; one for each unique scale.
% hashtable: a hashtable for the filters and ranges tables. The hash
% function for scale is round(scale,2)*100.
%
%% ----------------------------------------------------------------------------
kNUMPOINTS = int32(size(points,1));

%% Select orientation for each point..
% and describe with a row vector of length 2 [theta, phi].
%[orientation,isvalid] = assignorientation(V,points,scale,isvalid);
orientation(kNUMPOINTS,2) = single(0); % U-SURFs

%% Grab oriented chunks from the Volume.

% Determine the radius of the inital chunk accounting for rotation.
regionsize = computeregionsize(windowsize,scale);
kRADIUS = regionradius(windowsize,scale,orientation);
[x,y,z] = size(V);

region = cell(kNUMPOINTS,1);
parfor i = 1:kNUMPOINTS
    if isvalid(i)
    % See if the rotated region intersects the edge of the volume.
    if(sum(1 > (points(i,:) - kRADIUS(i,:)),2) > 0 || sum([x,y,z] < (points(i,:) + kRADIUS(i,:)),2) > 0)
        %fprintf(1,'\nREGION %i OUT OF BOUNDS.', i);
        isvalid(i) = false;
    else
        % Calculate the corners of the rotated region in V space.
        top = points(i,:) + kRADIUS(i,:);
        bot = points(i,:) - kRADIUS(i,:);
        
        % Cut out the chunk to be rotated.
        region{i} = V(bot(1):top(1),bot(2):top(2),bot(3):top(3));
        
        % Define the tranformation.
        T = maketform('affine',makeRmatrix(orientation));
        R = makeresampler('linear', 'fill');
        region{i} = tformarray(region{i},T,R,[1,2,3],[1,2,3],...
                      double([regionsize(i),regionsize(i), regionsize(i)]),[],NaN);
        
        % Debugging radius size.          
        assert(sum(sum(sum(isnan(region{i})))) == 0);
    end
    end
end
clear orientation kRADIUS;

%% Make a hashtable to look up filters and sampling spaces. -------------------

% Find unique scales and sort smallest to largest.
% Don't hash scales that have no valid points.
[uniquevalues,hashtable] = makehashtable(scale, isvalid); 
kTABLESIZE = length(uniquevalues);
regionsize = computeregionsize(windowsize,uniquevalues);
subregionsize = computesubregionsize(windowsize,uniquevalues);

% Create lookup tables for filters and sampling spaces.
ranges(kTABLESIZE,5) = int16(0);
filters = cell(kTABLESIZE,1);
for i = 1:kTABLESIZE
    % Make a gaussian filter for this scale.
    filters{i} = gaussian3D(regionsize(i), 3.3*uniquevalues(i));
    
    % Decide where to sample each subregion for this scale.
    range = zeros(1,5,'uint8'); % integer operations
    range(3) = subregionsize(i)/2 + 1;
    range(5) = subregionsize(i) - int16(uniquevalues(i));
    range(1) = 1 + (uniquevalues(i));
    range(2) = mean([range(1),range(3)]);
    range(4) = mean([range(3),range(5)]);
    ranges(i,:) = range;
end
end
%% Helper Functions -----------------------------------------------------------

function [R] = makeRmatrix(orientation)

%% U-SURFS: All orientations are upright.
Rtheta = [   1,  0,  0,  0;
             0,  1,  0,  0;
             0,  0,  1,  0;
             0,  0,  0,  1   ];
        
Rphi =   [   1,  0,  0,  0;
             0,  1,  0,  0;
             0,  0,  1,  0;
             0,  0,  0,  1   ];
%% Non-USURFS
             
R = Rtheta*Rphi;
        
end

function [radius] = regionradius(windowsize,scale,orientation)
    radius = int16(ceil(1.41*windowsize.*scale/2));
    radius = int16(ceil(windowsize.*scale/2));
end

function [rsize] = computeregionsize(windowsize,scale)
    rsize = int16(windowsize.*scale)+1;
end

function [srsize] = computesubregionsize(windowsize,scale)
    srsize = int16(ceil((windowsize.*scale+1)/3));
end

function [uniquevalues,hashtable] = makehashtable(scale, isvalid)
%MAKEHASHTABLE makes a hashtable for the valid values of scale using the
% hashfunction h = scale*100;
%
% INPUT
% scale: scale values for the points rounded to 2 decimal places.
% isvalid: a boolean array of whether each point had a valid descriptor.
%
% OUTPUT
% hashtable: a table where element 124 contains the index of scale 1.24 in
% the array UNIQUEVALUES.
% uniquevalues: a list of the unique scales in SCALE.
%
%% -----------------------------------------------------------------------
uniquevalues = unique(scale.*isvalid);
if(uniquevalues(1) == 0), uniquevalues = uniquevalues(2:end); end;
hashedvalues = int16(uniquevalues.*100);
kNUMUNIQUEVALUES = length(uniquevalues);

if(kNUMUNIQUEVALUES >= 256)
    error('Hashtable bitdepth not suffient for number of unique scales.');
end

hashtable(hashedvalues(kNUMUNIQUEVALUES),1) = uint8(0);
for i = 1:kNUMUNIQUEVALUES
    hashtable(hashedvalues(i)) = i;
end
end
