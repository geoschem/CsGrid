function [ midLL, midXYZ ] = cellCenter( varargin )
%CELLCENTER Get the lat/lon center of a cell
%   Expects 4 lat-lon pairs as inputs (cell corners, radians)

if length(varargin) == 1
    assert(all(size(varargin{1}) == [2,4]),'cellCenter:badArg1','Need 4 arguments (cell corners)');
    LLCorner = varargin{1};
else
    assert(length(varargin)==4,'cellCenter:badArg2','Need 4 arguments (cell corners)');
    LLCorner = zeros(2,4);
    for iPt = 1:4
        LLCorner(:,iPt) = varargin{iPt};
    end
end
midXYZ = zeros(3,1);
for iPt = 1:4
    midXYZ = midXYZ + LL2XYZ2(LLCorner(1,iPt),LLCorner(2,iPt));
end

% Normalize
sumVal = sum(midXYZ.^2.0);
if sumVal > 0
    midXYZ = midXYZ ./ sqrt(sumVal);
end

% Convert back to lat-lon
midLL = zeros(1,2);
[midLL(1),midLL(2)] = XYZ2LL2(midXYZ);

end

