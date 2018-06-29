function [ res ] = logRand( a, b, n)
%LOGRAND Summary of this function goes here
%   Detailed explanation goes here
    LA = log10(a); 
    LB = log10(b);
    res = 10.^(LA + (LB-LA) * rand(1,n));
end

