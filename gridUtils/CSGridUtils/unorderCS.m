function [ CSData ] = unorderCS( CSData, isGCHP )
%UNORDERCS Reverses the "reorderCS" operation
%   Detailed explanation goes here

% Assume data is order GCHP-style
if nargin < 2
    isGCHP = true;
end
nPerSide = size(CSData,1);
if isGCHP
    faceOrder = [1,2,6,4,5,3];
    faceRot = [0,0,90,180,180,90];
else
    faceOrder = 1:6;
    faceRot = zeros(1,6);
end

CS3D = CSData;
CSData = zeros(nPerSide,6*nPerSide);
for iFace = 1:6
    temp2D = CS3D(:,:,iFace);
    % Change WRITE order
    readFace = faceOrder(iFace);
    faceLo = (nPerSide * (readFace-1)) + 1;
    faceHi = faceLo + nPerSide - 1;
    currRot = faceRot(iFace);
    n90 = round(currRot/90);
    if n90 ~= 0
        temp2D = rot90(temp2D,-1 * n90);
    end
    CSData(:,faceLo:faceHi) = temp2D;
end

end

