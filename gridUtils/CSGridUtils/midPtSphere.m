function [ C, Cxyz ] = midPtSphere( A, B )
%MIDPTSPHERE Gives the mid-point (C) of 2 points (A and B) on a sphere
%   Arguments and output are [Lon,Lat]

A = LL2XYZ2(A(1),A(2));
B = LL2XYZ2(B(1),B(2));
Cxyz = midPt2(A,B);
C = zeros(1,2);
[C(1),C(2)] = XYZ2LL2(Cxyz);

end

function [v12] = midPt2(xyz1,xyz2)
%MIDPT2 Cartesian mid-point of 2 other points

v12 = xyz1 + xyz2;
v12 = v12./sqrt(sum(v12.^2));

end