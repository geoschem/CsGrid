function [ uOut, vOut ] = rotateLLCS( uIn, vIn, leftLL, baseLL, rightLL, topLL, centerLL, LLInput )
%ROTATELLCS Rotates a vector between lat-lon and cubed-sphere
%   Expect 
%   Inputs:
%       uIn:        Magnitude of u-component (aligned with source cell)
%       vIn:        Magnitude of v-component (aligned with source cell)
%       leftLL:     Mid-point of left edge of CS cell
%       baseLL:     Mid-point of bottom edge of CS cell
%       rightLL:    Mid-point of right edge of CS cell
%       topLL:      Mid-point of top edge of CS cell
%       centerLL:   Point at which the wind is defined
%       LLInput:    If true, rotate LL => CS. Otherwise, CS => LL.

nDimIn = 2;
allVec = {leftLL,baseLL,rightLL,topLL,centerLL};
for iVec = 1:5
    currVec = allVec{iVec};
    assert(numel(currVec) == nDimIn,'rotateLLCS:badBaseVec','Input vectors must have %i elements for this conversion',nDimIn);
    assert(max(abs(currVec)) <= 3*pi,'rotateLLCS:onlyRads','Input must be in radians');
end

% Calculate unit vector pointing from left to right
eLR = getUnitVecLL(rightLL,centerLL,leftLL);
% Do the same for bottom to top
eBU = getUnitVecLL(topLL,centerLL,baseLL);

eLon = zeros(3,1);
eLat = zeros(3,1);

eLon(1) =                   -sin(centerLL(1) - pi);
eLon(2) =                    cos(centerLL(1) - pi);
eLon(3) =  0.0                                    ;
eLat(1) = -sin(centerLL(2))* cos(centerLL(1) - pi);
eLat(2) = -sin(centerLL(2))* sin(centerLL(1) - pi);
eLat(3) =  cos(centerLL(2))                       ;

% or...
eLat = cross(LL2XYZ2(centerLL(1),centerLL(2)),eLon);

g11 = dot(eLR,eLon);
g12 = dot(eLR,eLat);
g21 = dot(eBU,eLon);
g22 = dot(eBU,eLat);

if LLInput
    uOut = uIn * g11 + vIn * g12;
    vOut = uIn * g21 + vIn * g22;
else
    uOut = ( uIn * g22 - vIn * g12)/(g11*g22 - g21*g12);
    vOut = (-uIn * g21 + vIn * g11)/(g11*g22 - g21*g12);
end

% Question: why is (uIn.*eLon + vIn.*eLat) != (uOut.*eLR + vOut.*eBU)??
% Answer:   eLon and eLat are not coplanar with eLR and eBU?!

end

function [eVec] = getUnitVecLL(startLL,midLL,endLL)
%GETUNITVECLL Get cartesian unit vector (on sphere) from 3 lat-lon points

endPt = LL2XYZ2(endLL(1),endLL(2));
midPt = LL2XYZ2(midLL(1),midLL(2));
startPt = LL2XYZ2(startLL(1),startLL(2));
eVec = endPt - startPt;

% Project onto sphere
dProd = dot(midPt,eVec);
eVec = eVec - (dProd.*midPt);

% Normalize
eVec = eVec./sqrt(sum(eVec.^2));

end