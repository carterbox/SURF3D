function [] = surfhessian3D( filter_size )
% SURFHESSIAN3D approximates the hessian matrix by using box filters as
% described in the SURF paper by H Bay, A Ess, T Tuytelaars, and L Van Gool
% in 2008.

% The filter size must be multiple of 3 and greater than 9 in size.
assert( mod(filter_size, 3) == 0 && filter_size > 9 );

% The Hessian matrix elements in this order: D11, D22, D33, D12, D13, D23.
H = zeros(6,1);





end

%% TODO: Fix these second partial derivatives so that they taper like the 3D gaussians do. I.E. make the z direction non constant.

function [cube] = makeparital11(filter_size)
% Makes the 2nd partial derivative box filter with respect to direction 1
% and then direction 1. FILTER_SIZE is a scalar because the box is square.

% Determine the characterisitc length.
l0 = filter_size/3;

% Make the stripes that define the single direction second derivative.
top = ones(l0,l0*2 - 1);
middle = -2.*ones(l0,l0*2 - 1);

% Assemble the stripes and pad into a cube.
box = cat(1, top, middle, top);
box = padarray(box, [0,0,l0 - 1], 'replicate'); 
cube = padarray(box, [0,(l0 + 1)/2,(l0 + 1)/2]);

% Make sure the box is square.
assert(size(cube,1) == filter_size && size(cube,2) == filter_size,...
    size(cube,3) == filter_size);
end

function [cube] = makepartial12(filter_size)
% Makes the 2nd partial derivative box filter with respect to direction 1
% and then direction 2.

% Determine the characterisitc length.
l0 = filter_size/3;

% Make the patches that define the filter.
top = cat(2, zeros(l0,(l0 - 1)/2), ones(l0), zeros(l0,1),...
    -ones(l0), zeros(l0,(l0 - 1)/2));
bottom = -top;

% Assemble the filter
box = cat(1, zeros((l0 - 1)/2, filter_size), top, zeros(1,filter_size),...
    bottom, zeros((l0 - 1)/2));
cube = padarray(box, [0,0,l0], 'replicate');
cube = padarray(box, [0,0,(l0 - 1)/2]);

% Make sure the box is square.
assert(size(cube,1) == filter_size && size(cube,2) == filter_size,...
    size(cube,3) == filter_size);
end