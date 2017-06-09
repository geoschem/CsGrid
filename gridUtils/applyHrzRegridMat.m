function [ outData ] = applyHrzRegridMat( inData, xDataObj )
%APPLYHRZREGRIDMAT Applies a regrid object to input data
%   Regridding object must be generated using calcHrzRegridMat. If the data
%   is sized like the output rather than input of the regrid object, the
%   inverse transform will be applied.

% Retrieve the relevant data
xRegrid = xDataObj.xRegrid;
gridOut = xDataObj.gridOut;

% Figure out the size of the input grid
inSize = size(inData);

% Make sure the shapes conform
nElData = inSize(1)*inSize(2);
nElOut = size(xRegrid,2);

if nElData ~= size(xRegrid,1)
    % Do we want to use the inverse transform?
    if nElData == nElOut
        % Invert transform
        xRegrid = xRegrid';
        gridOut = xDataObj.gridIn;
    else
        error('applyHrzRegridMatGeneric:incompatibleTransform',...
            'Input data is incorrectly sized for the given transform to be run either forwards or backwards');
    end
end

% Calculate weighting
testVec = ones(1,nElData);
weightVec = testVec*xRegrid;

outSize = inSize;
outSize(1:2) = gridOut;

% Sparse data multiplication requires double precision
reCastIn = ~isa(inData,'double');

if nElData == prod(inSize)
    % Only 1 slice - much more efficient to hardcode this case
    % Reshape input data into a vector
    inVec = reshape(inData,1,[]);
    if reCastIn
        inVec = cast(inVec,'double');
    end
    % Calculate as a vector first
    outVec = inVec*xRegrid;
    % Normalize
    outData = reshape(outVec./weightVec,gridOut);
else
    % Can handle arrays of N > 2 dimensions for any N
    outData = zeros(outSize,'like',inData);
    nSlice = prod(inSize(3:end));
    for iSlice = 1:nSlice
        % Reshape input data into a vector
        inVec = reshape(inData(:,:,iSlice),1,[]);
        if reCastIn
            inVec = cast(inVec,'double');
        end
        % Calculate as a vector first
        outVec = inVec*xRegrid;
        % Normalize
        outData(:,:,iSlice) = reshape(outVec./weightVec,gridOut);
    end
end

clear weightVec;

end

