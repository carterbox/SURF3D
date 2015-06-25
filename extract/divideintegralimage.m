function [subregion, subfilter] = divideintegralimage(J,filter,split,kNUMSUBREGIONS)
%DIVIDEINTEGRALIMAGE divides each region and its accompanying filter into
% subregions.
%
% INPUTS
% J (double): the integral image
% filter: the filter used to weight the haar wavelet responses.
% split: 1x3 vector telling how to split the volume in each direction.
%
% OUTPUTS
%
%% -----------------------------------------------------------------------
assert(size(J,1)-1 == size(filter,1) && size(J,2)-1 == size(filter,2) && size(J,3)-1 == size(filter,3)); 

% Divide the subfilter into n ranges along each direction.
[x,~,~] = size(J);
xrange = overlappingdivide(1,x,split(1));
yrange = xrange;%overlappingdivide(1,y,split(2));
zrange = xrange;%overlappingdivide(1,z,split(3));

% Add one to the ranges for the integral image because it's one bigger.
Jxrange = cat(2, xrange, xrange(:,end)+1);
Jyrange = Jxrange;%cat(2, yrange, yrange(:,end)+1);
Jzrange = Jxrange;%cat(2, zrange, zrange(:,end)+1);

% Generate a grid for the coordinates of each subregion.
[X,Y,Z] = ndgrid(1:split(1),1:split(2),1:split(3));

subregion = cell(kNUMSUBREGIONS,1);
subfilter = cell(kNUMSUBREGIONS,1);
parfor i = 1:kNUMSUBREGIONS
    subfilter{i} = filter(xrange(X(i),:),yrange(Y(i),:),zrange(Z(i),:));
    subregion{i} = J(Jxrange(X(i),:),Jyrange(Y(i),:),Jzrange(Z(i),:));
end

end

function ranges = overlappingdivide(lo,hi,n)
%OVERLAPPINGDIVIDE divides the range [lo,hi) into N equal and sometimes
% overlapping regions.
% 
% INPUTS
% lo (double):
% hi (double):
% n (int):
%
%% -----------------------------------------------------------------------
if(n < 1), error('Cannot divide into less than one group'); end;
if(hi-lo < n), error('Range must be greater or equal to %i.', n); end;

width = ceil((hi-lo)/double(n)); % The size of each subrange.

ranges(n,width) = uint8(0);
i = uint8(1); j = uint8(n); placeleft = true;
left = lo; right = hi;
while lo < hi
    if placeleft
        % Place a subrange on the left side
        left = lo+width;
        % Shift to center if overlaping
        if(left > right)
            lo = lo - (left-right)/2;
            left = lo+width;
        end
        ranges(i,:) = lo:left-1;
        lo = left;
        i = i+1;
    else
        % Place a subrange on the right size.
        right = hi-width;
        % Shift to center if overlaping
        if(left > right)
            hi = hi + (left-right)/2;
            right = hi-width;
        end
        ranges(j,:) = right:hi-1;
        hi = right;
        j = j-1;
    end
    placeleft = ~placeleft;
end
end