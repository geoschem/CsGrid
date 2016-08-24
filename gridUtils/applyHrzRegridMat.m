function [ dataOut ] = applyHrzRegridMat( dataIn, xMatLon, xMatLat )
%APPLYHRZREGRIDMAT Applies horizontal regridding transformations

% Get data sizes
nLatIn = size(xMatLat,2);
nLatOut = size(xMatLat,1);
nLonOut = size(xMatLon,1);
dataMid = zeros(nLonOut,nLatIn,'like',dataIn);

if issparse(xMatLon)
    [outEntry,inEntry,xVal] = find(xMatLon);
    for iEntry = 1:length(xVal)
        iOut = outEntry(iEntry);
        % NOTE: Longitudes can loop around!
        iIn = mod(inEntry(iEntry)-1,size(dataIn,1))+1;
        dataMid(iOut,:) = dataMid(iOut,:) + xVal(iEntry)*dataIn(iIn,:);
    end
else
    for iLat = 1:nLatIn
        % Regrid horizontally
        dataMid(:,iLat) = xMatLon*dataIn(:,iLat);
    end
end
clear dataIn;
dataOut = zeros(nLonOut,nLatOut,'like',dataMid);
if issparse(xMatLat)
    [outEntry,inEntry,xVal] = find(xMatLat);
    for iEntry = 1:length(xVal)
        iOut = outEntry(iEntry);
        iIn = inEntry(iEntry);
        dataOut(:,iOut) = dataOut(:,iOut) + xVal(iEntry)*dataMid(:,iIn);
    end
else
    for iLon = 1:nLonOut
        dataOut(iLon,:) = xMatLat*dataMid(iLon,:)';
    end
end

end

