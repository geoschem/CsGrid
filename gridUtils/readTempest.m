function [ xData ] = readTempest( TFPath )
% READTEMPEST Read Tempest NetCDF tile file
%             regridding information to generate an xData struct array 
%             that can generate a GCHP-compatible tile file when passed
%             to writeTileFile().
%
%    Input:  (1) Tfile: Tempest output netcdf file (full path)
%
%    Output: (1) xData: struct array with same format created when reading
%                       a GMAO tile file using readTileFile.m
%
% NOTES: 
%    (1) There is currently a stretching discrepancy between cubed sphere
%        data sets created with GMAO vs Tempest mapping information.
%    (2) Tempest output currently does not include half-polar options,
%        so tilefiles created from xData will be of form DExPE_CFx6C.
%
% Lizzie Lundgren, 10/7/16

%------------------------
% Check that file exists
%-----------------------
assert(fastExist(TFPath),'readTempest:fileNotFound', ...
       'Tempest file does not exist');

% Get the number of elements in each array
TFInfo = ncinfo(TFPath);
nElIn = ncread(TFPath,'src_grid_dims');
nElOut = ncread(TFPath,'dst_grid_dims');

% Lat-lon will be 2 elements (nLon x nLat), cubed-sphere will be 1 element
% (N*N*6)
isCS = false(2,1);
[isCS(1),gDimsIn,gNameIn] = parseDims(nElIn,TFPath,'a');
[isCS(2),gDimsOut,gNameOut] = parseDims(nElOut,TFPath,'b');

% Read relevant data
colData = ncread(TFPath,'col');
rowData = ncread(TFPath,'row');
SData = ncread(TFPath,'S');

% Create xData
xData = struct( 'Name',{gNameIn;gNameOut},...
                'NX',{gDimsIn(1);gDimsOut(1)},...
                'NY',{gDimsIn(2);gDimsOut(2)},...
                'II',cell(2,1),...
                'JJ',cell(2,1),...
                'W',SData );

% The ind2sub function is much faster when called on a vector
%{
% How many points do we have?
numPoints = length(colData);
for i=1:numPoints
    % Each data point links:
    %   Input cell [x(1).II(i),x(1).JJ(i)]
    % to
    %   Output cell [x(2).II(i),x(2).JJ(i)]
    % with a factor
    %   Weight W(i)
    [xData(1).II(i), xData(1).JJ(i)] = ind2sub(gDimsIn, colData(i));
    [xData(2).II(i), xData(2).JJ(i)] = ind2sub(gDimsOut, rowData(i)); 
end
%}
[xData(1).II,xData(1).JJ] = ind2sub(gDimsIn, colData);
[xData(2).II,xData(2).JJ] = ind2sub(gDimsOut,rowData);

% Change the cubed sphere indexing conventions to match GMAO
faceRemapping =  [4 5 1 2 6 3];
for iGrid = 1:2
    if isCS(iGrid)
        % Swap faces
        xData(iGrid) = swapCSFaces( xData(iGrid), faceRemapping);

        % Flip face 6 in both dimensions
        xData(iGrid) = flipCSFaces( xData(iGrid), 6, '2d');

        % Transpose and flip (II only) faces 3-5
        xData(iGrid) = transposeCSFaces( xData(iGrid), [3 4 5]);
        xData(iGrid) = flipCSFaces( xData(iGrid), [3 4 5], 'II');
    end
end
end

function [isCS,gDims,gName] = parseDims(nEl,fPath,gridChar)
isCS = numel(nEl) == 1;
if isCS
    nCS = sqrt(cast(nEl/6,'double'));
    assert(abs(nCS - round(nCS)) < 1e-10,'readTempestGeneric:badCSGrid',...
        'CS grid has an element count which cannot be parsed');
    gDims = [nCS,nCS*6];
    gName = sprintf('CF%04ix6C',nCS);
else
    gDims = nEl;
    % Need to figure out quite a lot of information. Retrieve the vertices
    % describing the grid for whichever is the lat-lon. Each cell is
    % described by 4 vertices, and the array has size [4,nCells]
    % SDE 2017-06-10: The routine below is approximate, but very fast. A
    % more rigorous alternative routine would be to read all vertex data
    % from the file; however, this takes a very long time!
    %latBnds = ncread(fPath,['yv_',gridChar]);
    latCell = ncread(fPath,['yv_',gridChar],[1,1],[Inf,6]);
    lonCell = ncread(fPath,['xv_',gridChar],[1,1],[Inf,6]);
    dLon = median(max(lonCell,[],1)-min(lonCell,[],1));
    % Assume constant spacing in longitude
    % Are we on a subgrid?
    if nEl(1) * dLon < 359
        DStr = 'UU';
        PStr = 'UU';
    else
        % Currently the algorithm below only works some of the time.
        % A workaround is now included that checks the PE/PC and DE/DC strings
        % against the tempest filename which should contain the correct ones.
        % If the positioning is inconsistent then the value is over-written 
        % with what is in the filename (ewl, 12/15/17)
        [filepath,name,ext] = fileparts(fPath);
        msg1 = ['*** WARNING ***\n%s positioning is not correctly ' ...
               'interpreted in readTempest.m.\nResetting to ' ...
               'positioning in tempest filename...\n'];
        msg2 = ['*** WARNING ***\nPositioning is not in tempest file' ...
                ' to validate interpretation by readTempest.m.\n' ...
                'Beware that sometimes PE vs PC and DE vs DC are wrong.\n'];
	PEinName = contains(name,'PE');
	PCinName = contains(name,'PC');
	DEinName = contains(name,'DE');
	DCinName = contains(name,'DC');
        if ~PEinName && ~PCinName fprintf(msg2), end
        if ~DEinName && ~DCinName fprintf(msg2), end

        % If the first latitude delta is half the size of the next 4, we have a
        % "pole-centered" grid
        %dLat = latBnds(2,:) - latBnds(1,:);
        % Take the delta across the vertices
        %dLat = max(abs(diff(latBnds,[],1)),[],1);
        dLat = max(abs(diff(latCell,[],1)),[],1);
        dLatP = dLat(1);
        if length(dLat) > 5
            stopPt = 4;
        else
            stopPt = length(dLat) - 1;
        end
        dLatNP = mean(dLat(2:stopPt));
        if abs(dLatNP - (2*dLatP)) < 1e-6
            PStr = 'PC';
            if PEinName fprintf(msg1,'Latitude'), PStr = 'PE';, end
        else
            PStr = 'PE';
            if PCinName fprintf(msg1,'Latitude'), PStr = 'PC';, end
        end

        % Dateline-centered or pole-centered?
        %lonBnds = ncread(fPath,'lon_bnds');
        %lonBnds = reshape(ncread(fPath,['xv_',gridChar]),[],1);
        lonBnds = mod(lonCell(1) + cast(1:(nEl(1)+1),'double').*dLon,360);
        % If any of the longitude bounds are ~ 0, we are dateline-edged
        lonMis = min(lonBnds);
        if lonMis < 1e-8
            DStr = 'DE';
            if DCinName fprintf(msg1,'Longitude'), DStr = 'DC';, end
        else
            DStr = 'DC';
            if DEinName fprintf(msg1,'Longitude'), DStr = 'DE';, end
        end
    end
    gName = sprintf('%s%04ix%s%04i',DStr,nEl(1),PStr,nEl(2));
end
end