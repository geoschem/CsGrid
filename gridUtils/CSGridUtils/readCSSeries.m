function [ CSData ] = readCSSeries( targDir, varName, colName, varargin )
%READCSSERIES Read in a timeseries of data from NetCDF files

if ischar(varName)
    varName = {varName};
end
if ischar(colName)
    colName = {colName};
end
assert(iscell(varName) && iscell(colName),'readCSSeries:nonCellArg','Input arguments must be cell arrays');
nVar = length(varName);
if nVar > length(colName)
    if length(colName) == 1
        colName = repmat(colName,[nVar,1]);
    else
        error('readCSSeries:requestMismatch','Each variable requested must have a collection associated with it');
    end
end

% Check that the input directory exists!
assert(exist(targDir,'dir') > 0,'readCSSeries:missingDir','Input directory ''%s'' not found');

inParse = inputParser();
addParameter(inParse,'dateRange',[-Inf,+Inf],@(x)validateattributes(x,{'numeric'},{'numel',2}));
addParameter(inParse,'levRange',[],@(x)validateattributes(x,{'numeric'},{'positive','integer'}));
parse(inParse,varargin{:});

% First figure out common data: number of samples and grid sizes
unqCol = unique(colName);
nCollection = length(unqCol);
nFiles = zeros(1,nCollection);
fileNames = cell(1,nCollection);

% If a specific date range was requested, remove entries not in
% that range
dRange = inParse.Results.dateRange;
dMin = dRange(1);
dMax = dRange(2);

% Subrange of levels?
levRange = inParse.Results.levRange;
if isempty(levRange)
    nLev = 0;
else
    if isscalar(levRange)
        minLev = levRange;
        maxLev = levRange;
    else
        maxLev = max(levRange);
        minLev = min(levRange);
        assert(length(levRange) == (maxLev+1-minLev),...
            'readCSSeries:badLevRange','Level range must be continuous or bounding');
    end
    nLev = maxLev + 1 - minLev;
end

% Get the time samples available in each file
% This allows for uneven spacing, and for different timing between files
tCell = cell(nCollection,2);
outLength = zeros(nCollection,1);
for iCollection = 1:nCollection
    fullList = dir(fullfile(targDir,sprintf('*%s*.nc4',colName{iCollection})));
    fullList = {fullList.name};
    nFiles(iCollection) = length(fullList);
    fileNames{iCollection} = fullList;
    % Figure out the time vector in each
    tCellTemp = cell(nFiles(iCollection));
    nSamples = 0;
    nKeep = 0;
    for iFile = 1:nFiles(iCollection)
        currFile = fullfile(targDir,fileNames{iCollection}{iFile});
        dStart = ncreadatt(currFile,'time','begin_date'); % yyyymmdd
        tStart = ncreadatt(currFile,'time','begin_time'); % hhmmss
        tStart = datenum(sprintf('%08i %06i',dStart,tStart),'yyyymmdd HHMMSS');
        tVecTemp = cast(ncread(currFile,'time'),'double');
        % Assume it's in minutes
        tCellTemp{iFile} = tStart + (tVecTemp./(60.0*24.0));
        nSamples = nSamples + length(tVecTemp);
        nKeep = nKeep + sum(tVecTemp >= dMin & tVecTemp <= dMax);
    end
    outLength(iCollection) = nKeep;
    tVecFull = zeros(nSamples,1);
    tVecSrc = zeros(nSamples,1,'int32');
    hiPt = 0;
    for iFile = 1:nFiles(iCollection)
        tVecTemp = tCellTemp{iFile};
        loPt = hiPt + 1;
        hiPt = hiPt + length(tVecTemp);
        tVecFull(loPt:hiPt) = tVecTemp;
        tVecSrc(loPt:hiPt) = iFile;
    end
    [tVecFull,reorderVec] = sort(tVecFull);
    tVecSrc = tVecSrc(reorderVec);
    tCell{iCollection,1} = tVecFull;
    tCell{iCollection,2} = tVecSrc;
end

CSData = struct();
for iVar = 1:nVar
    % Add the time data
    currCol = colName{iVar};
    currVar = varName{iVar};
    safeName = matlab.lang.makeValidName(currVar);
    iCollection = find(strcmpi(currCol,unqCol));
    tVec = tCell{iCollection,1};
    tSrc = tCell{iCollection,2};
    whichFile = unique(tSrc);
    nRead = length(whichFile);
    gridReady = false;
    hiPt = 0;
    nKeepCol = outLength(iCollection);
    for iRead = 1:nRead
        % Which file are we reading from?
        iFile = whichFile(iRead);
        srcFile = fileNames{iCollection}{iFile};
        currFile = fullfile(targDir,srcFile);
        % Which of this file's samples do we actually want?
        tVecFile = tVec(tSrc == iRead);
        tKeep = (tVecFile >= dMin) & (tVecFile <= dMax);
        % Do we want a subrange of levels?
        varInfo = ncinfo(currFile,currVar);
        dimNames = {varInfo.Dimensions.Name};
        dimLen = [varInfo.Dimensions.Length];
        iTime = find(strcmpi(dimNames,'time'));
        if nLev > 0
            iLev = find(strcmpi(dimNames,'lev'));
            if isempty(iLev)
                readOpts = {};
            else
                readOpts = cell(2,1);
                readOpts{1} = ones(size(dimLen));
                readOpts{1}(iLev) = minLev;
                readOpts{2} = inf(size(dimLen));
                readOpts{2}(iLev) = nLev;
            end
        else
            readOpts = {};
        end
        baseData = ncread(currFile,currVar,readOpts{:});
        % Is the output grid allocated?
        if ~gridReady
            newSize = size(baseData);
            newSize(end) = nKeepCol;
            CSData.(safeName).data = zeros(newSize);
            targLoc = cell(1,length(newSize));
            for iLoc = 1:length(newSize)
                targLoc{iLoc} = ':';
            end
            CSData.(safeName).tVec = zeros(nKeepCol,1);
            gridReady = true;
        end
        loPt = hiPt + 1;
        hiPt = hiPt + sum(tKeep);
        assert(hiPt<=nKeepCol,'readCSSeries:timeMismatch','Could not process time vectors');
        targLoc{iTime} = loPt:hiPt;
        CSData.(safeName).data(targLoc{:}) = baseData;
        CSData.(safeName).tVec(loPt:hiPt) = tVecFile(tKeep);
    end
end

end