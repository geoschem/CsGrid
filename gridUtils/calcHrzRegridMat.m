function [ xDataObj ] = calcHrzRegridMat( gSpecIn, gSpecOut, TFDir )
%CALCHRZREGRIDMAT Retrieves or calculates conservative regridding weights
% Can accept any combination of lat-lon grid specifications and
% cubed-sphere specifier strings ('C180'). For anything other than a
% conventional LL -> LL regrid, a tile file must be present in "TFDir.
% If "TFDir" is not given, it is assumed that the relevant data can be
% found in GridData/TileFiles.

% Validate inputs
inParse = inputParser;
inParse.addRequired('gSpecIn',@(x)isa(x,'gridSpec') || ischar(x) && strncmpi(x,'c',1));
inParse.addRequired('gSpecOut',@(x)isa(x,'gridSpec') || ischar(x) && strncmpi(x,'c',1));
inParse.parse(gSpecIn,gSpecOut);
[nXIn,nYIn,CSIn,DCPCIn,bNameIn,NCNameIn] = determineGrid(gSpecIn);
[nXOut,nYOut,CSOut,DCPCOut,bNameOut,NCNameOut] = determineGrid(gSpecOut);

% Set up the output object
xDataObj.gridIn = [nXIn,nYIn];
xDataObj.gridOut = [nXOut,nYOut];
% Placeholder
xDataObj.xRegrid = [];

% Do we need a tile file?
needTF = (CSOut || CSIn);
if needTF
    if nargin < 3 || isempty(TFDir)
        TFDir = fullfile('GridData','TileFiles');
    end
    assert(exist(TFDir,'dir') == 7,'calcHrzRegridMatGeneric:noTFDir',...
        'Could not find tile files directory ''%s''',TFDir);
    % Generate binary and NetCDF names in each direction
    bNameI2O = sprintf('%s_%s.bin',bNameIn,bNameOut);
    bNameO2I = sprintf('%s_%s.bin',bNameOut,bNameIn);
    % NetCDF is a little trickier, as there is an assumption of only one
    % lat-lon grid in the transform
    if CSOut && CSIn
        NCNameI2O = sprintf('%s-to-%s_MAP_GMAO.nc',NCNameIn,NCNameOut);
        NCNameO2I = sprintf('%s-to-%s_MAP_GMAO.nc',NCNameOut,NCNameIn);
    else
        if CSOut
            DCPCStr = DCPCIn;
        else
            DCPCStr = DCPCOut;
        end
        NCNameI2O = sprintf('%s-to-%s_MAP_%s_GMAO.nc',NCNameIn,NCNameOut,DCPCStr);
        NCNameO2I = sprintf('%s-to-%s_MAP_%s_GMAO.nc',NCNameOut,NCNameIn,DCPCStr);
    end
    % Add the root to everything
    bNameI2O = fullfile(TFDir,bNameI2O);
    bNameO2I = fullfile(TFDir,bNameO2I);
    NCNameI2O = fullfile(TFDir,NCNameI2O);
    NCNameO2I = fullfile(TFDir,NCNameO2I);
    % Does the I2O tile file exist?
    if (fastExist(NCNameI2O))
        % Use Tempest NetCDF data (best)
        xData = readTempest(NCNameI2O);
    elseif (fastExist(bNameI2O))
        xData = readTileFile(bNameI2O);
    else
        % Note: it doesn't actually matter that the regridding
        % object is oriented incorrectly. This is handled in the
        % applyHrzRegridMatGeneric routine
        if (fastExist(NCNameO2I))
            xData = readTempest(NCNameO2I);
        elseif (fastExist(bNameO2I))
            xData = readTileFile(bNameO2I);
        else
            % Tell the user what Tempest settings to use
            fprintf('%s%s%s\n',repmat('=',[1,20]),' ERROR ',repmat('=',[1,20]));
            fprintf('Tile file could not be found.\n');
            if xor(CSOut,CSIn)
                if CSOut
                    nCS = nXOut;
                    gSpecLL = gSpecIn;
                else
                    nCS = nXIn;
                    gSpecLL = gSpecOut;
                end
                nLon = gSpecLL.nLon;
                nLat = gSpecLL.nLat;
                isNested = gSpecLL.isNested;
                if isNested
                    DCArg = 'false';
                    PCArg = 'false';
                else
                    if gSpecLL.center180
                        DCArg = 'true';
                    else
                        DCArg = 'false';
                    end
                    if gSpecLL.halfPolar
                        PCArg = 'true';
                    else
                        PCArg = 'false';
                    end
                end
                fprintf('Relevant tile file can be generated as\nfollows using Tempest:\n');
                if CSIn
                    fprintf(' => isC2L:  true\n');
                    fprintf(' => isL2C:  false\n');
                elseif CSOut
                    fprintf(' => isC2L:  false\n');
                    fprintf(' => isL2C:  true\n');
                end
                fprintf(' => nLon:   %i\n',nLon);
                fprintf(' => nLat:   %i\n',nLat);
                fprintf(' => nC:     %i\n',nCS);
                fprintf(' => isDC:   %s\n',DCArg);
                fprintf(' => isPC:   %s\n',PCArg);
                fprintf(' => isGMAO: true\n');
                if isNested
                    lonStart = mod(gSpecLL.lonEdge(1),360);
                    lonStop = mod(gSpecLL.lonEdge(end),360);
                    latStart = gSpecLL.latEdge(1);
                    latStop = gSpecLL.latEdge(end);
                    fprintf('Regional mask settings are also required:\n');
                    fprintf(' => region='' --lon_begin %f --lon_end %f --lat_begin %f --lat_end %f\n',...
                        lonStart,lonStop,latStart,latStop);
                end
            else
                % Grid names are strings
                fprintf('Tempest must be run to generate a %s to %s regridding file\n',gSpecIn,gSpecOut);
            end
            fprintf('Output file should be named %s and moved to %s\n',...
                NCNameI2O,TFDir);
            fprintf('%s%s%s\n',repmat('=',[1,20]),' ERROR ',repmat('=',[1,20]));
            error('calcHrzRegridMatGeneric:TileFileMissing',...
                'Could not find any tile files for the given configuration');
        end
    end
    % Convert the tile file to a regridding map
    xDataObj.xRegrid = genRegridObj(xData);
else
    % Extract the important data
    lonEdgeIn = gSpecIn.lonEdge;
    latEdgeIn = gSpecIn.latEdge;
    lonEdgeOut = gSpecOut.lonEdge;
    latEdgeOut = gSpecOut.latEdge;
    
    inSub = max(lonEdgeIn) - min(lonEdgeIn) < 359;
    if inSub
        % Input grid is not global
        lonEdgeIn = [lonEdgeIn(end)-360, lonEdgeIn(:)'];
    end
    outSub = max(lonEdgeOut) - min(lonEdgeOut) < 359;
    if outSub
        % Output grid is not global
        lonEdgeOut = [lonEdgeOut(end)-360, lonEdgeOut(:)'];
    end
    
    % Call the regridding function to do the actual work
    [xLon,xLat] = calcLLXMat(lonEdgeIn,latEdgeIn,...
        lonEdgeOut,latEdgeOut);
    
    % Remove the padding cell if the range was not global
    if inSub
        xLon = xLon(:,2:end);
    end
    if outSub
        xLon = xLon(2:end,:);
    end

    % Start conversion to a generic regridding object
    % This is equivalent to genRegridObj for tile file data
    nElIn = nXIn*nYIn;
    nElOut = nXOut*nYOut;

    % Loop over all elements in the xLon/xLat matrices
    [xLonI,xLonJ,xLonW] = find(xLon);
    [xLatI,xLatJ,xLatW] = find(xLat);
    nLonX = length(xLonI);
    nLatX = length(xLatI);
    fromVec = zeros(nLonX*nLatX,1);
    toVec = zeros(nLonX*nLatX,1);
    WVec = zeros(nLonX*nLatX,1);
    iPt = 0;
    for iLonPt = 1:nLonX
        % Current longitude match
        lonTo = xLonI(iLonPt);
        lonFrom = xLonJ(iLonPt);
        lonWeight = xLonW(iLonPt);
        for iLatPt = 1:nLatX
            % Current latitude match
            latTo = xLatI(iLatPt);
            latFrom = xLatJ(iLatPt);
            latWeight = xLatW(iLatPt);
            % Add to indexing arrays
            iPt = iPt + 1;
            fromVec(iPt) = sub2ind([nXIn,nYIn],lonFrom,latFrom);
            toVec(iPt) = sub2ind([nXOut,nYOut],lonTo,latTo);
            WVec(iPt) = lonWeight*latWeight;
        end
    end
    % Output
    xDataObj.xRegrid = sparse(fromVec,toVec,WVec,nElIn,nElOut);
end

end

function [nX,nY,isCS,DCPCStr,bName,NCName] = determineGrid(gSpec)
if ischar(gSpec)
    % Cubed-sphere grid
    [nX,nY,bName,NCName] = parseCS(gSpec);
    isCS = true;
    % Other options do not apply
    DCPCStr = '';
else
    %  Lat-lon grid
    isCS = false;
    [nX,nY,isPC,isDC,isUU,bName,NCName] = parseLL(gSpec);
    if isUU || ~(isDC || isPC)
        % Convention is to specify DEPE for non-global grids,
        % but it's not strictly necessary
        DCPCStr = 'DEPE';
    else
        if isDC
            DChar = 'C';
        else
            DChar = 'E';
        end
        if isPC
            PChar = 'C';
        else
            PChar = 'E';
        end
        DCPCStr = sprintf('D%sP%s',DChar,PChar);
    end
end
end

function [nX,nY,bName,NCName] = parseCS(CSString)
%PARSECS Read CS resolution integer from astring
nCS = str2double(CSString(2:end));
nX = nCS;
nY = nCS * 6;
% Name for binary files
bName = sprintf('CF%04ix6C',nCS);
% Name for NetCDF
NCName = sprintf('c%i',nCS);
end

function [nX,nY,isPC,isDC,isUU,bName,NCName] = parseLL(gSpec)
nX = gSpec.nLon;
nY = gSpec.nLat;
% Special considerations
% Half-polar ('pole-centered') grid?
isPC = gSpec.halfPolar;
% Dateline-centered grid?
isDC = gSpec.center180;
% Name for binary files
isUU = gSpec.isNested;
if isUU
    bName = genMAPLLLGridName([nX,nY]);
else
    bName = genMAPLLLGridName([nX,nY],isDC,isPC);
end
% Name for NetCDF
NCName = sprintf('lon%i_lat%i',nX,nY);
end
