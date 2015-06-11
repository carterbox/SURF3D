function [ n ] = euclidnorm( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n = sqrt(sum(sum(sum(A.*A))));

end


