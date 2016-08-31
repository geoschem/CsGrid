function [ outputData, xData ] = regridConservative( inputData, xData )
%REGRIDCONSERVATIVE Regrid data using a MAPL tile file
%   xData can be generated using "readTileFile".

global CSGridDir

XYIn = size(inputData);
XYIn = XYIn(1:2);
if isstruct(xData)
    % Determine which transformation we are doing
    XY1 = [xData(1).NX,xData(1).NY];
    XY2 = [xData(2).NX,xData(2).NY];
    if all(XYIn == XY1)
        iIn = 1;
        iOut = 2;
    elseif all(XYIn == XY2)
        iIn = 2;
        iOut = 1;
    else
        error('regridConservative:invalidInput','Input incorrectly sized for given regrid data');
    end
else
    % Try to load a file, if necessary
    if isa(xData,'gridSpec')
        % Try to figure out what the input and output data are
        gSpec = xData;
        XYOut = [gSpec.nLon,gSpec.nLat];
        outGrid = genLLGridName(XYOut,gSpec.center180,gSpec.halfPolar);
    elseif ischar(xData) && strncmpi(xData,'c',1)
        % Cubed-sphere?
        nPerSide = str2double(xData(2:end));
        outGrid = sprintf('CF%04ix6C',nPerSide);
    else
        error('regridConservative:badXData','Transfer parameters not recognized');
    end
    % Is the input data cubed-sphere?
    if (XYIn(2) == 6*XYIn(1))
        inGrid = sprintf('CF%04ix6C',XYIn(1));
    else
        % Assume lat-lon
        % If an even number of latitudes, assume that we are pole-edged;
        % otherwise, assume pole-centered. For now, always assume
        % dateline-centered (GMAO standard)
        inGrid = genLLGridName(XYIn,true,mod(XYIn(2),2) == 1);
    end
    % Try to read the file
    xData = sprintf('%s_%s.bin',inGrid,outGrid);
    iIn = 1;
    iOut = 2;
    xFile = fullfile(CSGridDir,'GridData','TileFiles',xData);
    if ~fastExist(xFile)
        xData = sprintf('%s_%s.bin',outGrid,inGrid);
        xFile = fullfile(CSGridDir,'GridData','TileFiles',xData);
        iIn = 2;
        iOut = 1;
        assert(fastExist(xFile),'regridConservative:tileFileNotFound',...
            'Tile file does not exist for this configuration (%s)',xData);
    end
    xData = readTileFile(xFile);
end

% Allow for up to 3 additional dimensions
% NOTE: This is completely arbitrary, mostly included as a safeguard
% against doing something really inefficient. This number can be increased,
% but you will need to code for the additional dimensions
maxExtraDims = 2;
nDimRegrid = ndims(inputData) - 2;
assert(nDimRegrid<=maxExtraDims,'regridConservative:tooManyDims',...
    'Conservative regridding function currently limited to a maximum of %i additional dimensions',...
    maxExtraDims);

inSize = size(inputData);
outSize = inSize;
outSize(1) = xData(iOut).NX;
outSize(2) = xData(iOut).NY;
outputData = zeros(outSize,'like',inputData);
outXYSize = outSize(1:2);
nD4 = size(inputData,4);
nD3 = size(inputData,3);
nPoints = length(xData(iOut).II);

nanValue = nan(1);
weightSum = zeros(outXYSize,'like',inputData);
for N = 1:nPoints
    IIn = xData(iIn).II(N);
    JIn = xData(iIn).JJ(N);
    inValue = inputData(IIn,JIn,:,:); % 23.5%, 0.721 s
    inValue(isnan(inValue)) = nanValue;
    II = xData(iOut).II(N);
    JJ = xData(iOut).JJ(N);
    W  = xData(iOut).W(N);
    outputData(II,JJ,:,:) = outputData(II,JJ,:,:) + W.*inValue; % 43%, 1.32s
    weightSum(II,JJ) = weightSum(II,JJ) + W;
end
% Manually impose NaNs rather than relying on a div-by-zero error
noValue = weightSum==0;
weightSum(noValue) = 1;
for iD4 = 1:nD4
    for iD3 = 1:nD3
        outSlice = outputData(:,:,iD3,iD4)./weightSum;
        outSlice(noValue) = nan;
        outputData(:,:,iD3,iD4) = outSlice;
    end
end

end

function gName = genLLGridName(XYVec,center180,halfPolar)
if halfPolar
    poleChar = 'C';
else
    poleChar = 'E';
end
if center180
    dateChar = 'C';
else
    dateChar = 'E';
end
gName = sprintf('D%s%04ixP%s%04i',dateChar,XYVec(1),poleChar,XYVec(2));
end
