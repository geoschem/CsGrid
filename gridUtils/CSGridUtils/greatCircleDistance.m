function [ gcDist ] = greatCircleDistance( LL1, LL2, radius )
%GREATCIRCLEDISTANCE Calculate the great circle distance between two points
%   LL1 and LL2 are 2-element vectors ([longitude, latitude] in rad)
%   radius is assumed to be 1 if not given

% This formula works too, but may as well use MATLAB builtin
%{
betaVal = asin(sqrt(sin((LL1(2)-LL2(2))/2.0).^2.0 + cos(LL1(2))*cos(LL2(2))...
    *sin((LL1(1)-LL2(1))/2.0).^2.0)).*2.0;
gcDist = betaVal*radius;
%}
LL1 = LL1 .* 180./pi;
LL2 = LL2 .* 180./pi;
gcDist = radius * (pi/180) * distance(LL1(2),LL1(1),LL2(2),LL2(1));

end

