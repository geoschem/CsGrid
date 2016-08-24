function [ patchVec ] = updateCSPatch( patchVec,CSData, GCHPData )
%UPDATECSPATCH Update colors on CS grid
%   Takes the output from plotCSLayer and updates the patch data
%   Only valid for "edge" output

if GCHPData
    faceOrder = [1,2,6,4,5,3];
    faceRot = [0,0,90,180,180,90];
else
    faceOrder = 1:6;
    faceRot = zeros(1,6);
end
nPerSideSq = length(get(patchVec(1),'cdata'));
nPerSide = sqrt(nPerSideSq);
assert(nPerSide==floor(nPerSide),'updateCSPatch:badSize','Patch data incorrectly sized');
if ismatrix(CSData)
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
    clear CS2D;
end

% Update the faces
for iFace = 1:6
    dataVec = nan(nPerSide*nPerSide,1);
    faceData = CSData(:,:,iFace);
    iIndex = 0;
    for iLon = 1:nPerSide
        for iLat = 1:nPerSide
            iIndex = iIndex + 1;
            dataVec(iIndex) = faceData(iLon,iLat);
        end
    end
    set(patchVec(iFace),'cdata',dataVec);
end

end

