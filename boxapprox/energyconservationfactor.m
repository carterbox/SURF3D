function [w] = energyconservationfactor()

Lxx = euclidnorm(gaussian3D(9,1.2,11));
Dxx = euclidnorm(makefilter11(9));
Lxy = euclidnorm(gaussian3D(9,1.2,12));
Dxy = euclidnorm(makefilter12(9));

w = Lxy/Lxx * Dxx/Dxy;

end

function [box] = makefilter11(filter_size)
% Makes a Dxx filter of size 9.

base = cat(1, ones(3,5), -2*ones(3,5), ones(3,5));
base = cat(3,base,base,base,base,base);
box = padarray(base,[0,2,2]);

end

function [box] = makefilter12(filter_size)
% Makes a Dxy filter of size 9.

base = cat(2, ones(3,3,7), zeros(3,1,7), -ones(3,3,7));
inside = cat(1, base, zeros(1,7,7), -base); 
box = padarray(inside,[1,1,1]);

end

function [ n ] = euclidnorm( A )
%EUCLIDNORM calculates the euclidian norm of an N-d matrix.
%   We can't use the Matlab built in because it doesn't handle N-d
%   matricies.

A2 = A.*A;
n = sqrt(sum(A2(:)));

end
