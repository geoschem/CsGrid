function [ CSData ] = reorderCS( CSData, isGCHP )
%REORDERCS Converts 2D (Nx6N) CS data to 3D (NxNx6)

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

CS2D = CSData;
CSData = zeros(nPerSide,nPerSide,6);
for iFace = 1:6
    % Change READ order
    readFace = faceOrder(iFace);
    faceLo = (nPerSide * (readFace-1)) + 1;
    faceHi = faceLo + nPerSide - 1;
    temp2D = CS2D(:,faceLo:faceHi);
    currRot = faceRot(iFace);
    n90 = round(currRot/90);
    if n90 ~= 0
        temp2D = rot90(temp2D,n90);
    end
    CSData(:,:,iFace) = temp2D;
end

end

