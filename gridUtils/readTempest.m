function [ xData ] = readTempest( Tfile, Nlon, Nlat, NX, dateline, polar )
% READTEMPEST Read Tempest netcdf file (lat/lon to cubed sphere) 
%             regridding information to generate an xData struct array 
%             that can generate a GCHP-compatible tile file when passed
%             to writeTileFile().
%
%    Input:  (1) Tfile: Tempest output netcdf file (full path)
%            (2) Nlon:  number longitudes
%            (3) Nlat:  number latitudes
%            (4) NX:    cubed sphere side length
%            (5) dateline: DE or DC
%            (6) polar: PE or PC
%
%    Output: (1) xData: struct array with same format created when reading
%                       a GMAO tile file using readTileFile.m
%
% NOTES: 
%    (1) There is currently a stretching discrepancy between cubed sphere
%        data sets created with GMAO vs Tempest mapping information.
%    (2) Tempest output currently does not include half-polar options,
%        so tilefiles created from xData will be of form 
%        {dateline}x{polar}_CFx6C.
%
% Lizzie Lundgren, 10/7/16

%------------------------
% Check that file exists
%-----------------------
assert(fastExist(Tfile),'readTempest:fileNotFound', ...
       'Tempest file does not exist');
fID = fopen(Tfile,'rb');
cleanupFn = onCleanup(@()(fclose(fID)));

%---------------------------------
% Read file (lat/lon -> cs)
%---------------------------------
L2C_L_1D = ncread( Tfile, 'col');
L2C_C_1D = ncread( Tfile, 'row');
L2C_S    = ncread( Tfile, 'S'  );

%---------------------------------
% Create xData
%---------------------------------
numPoints = length(L2C_L_1D);
lonStr = [repmat('0', 1, 4-length(num2str(Nlon))), num2str(Nlon)];
latStr = [repmat('0', 1, 4-length(num2str(Nlat))), num2str(Nlat)];
cStr   = [repmat('0', 1, 4-length(num2str(NX))),   num2str(NX)];
LLname = [dateline lonStr 'x' polar latStr];
CSname = ['CF' cStr 'x6C'];
xData = struct('Name', {LLname; CSname}, ...
	       'NX',   {Nlon; NX},       ...
	       'NY',   {Nlat; NX*6},     ...
	       'II',   cell(2,1),        ...
	       'JJ',   cell(2,1),        ...
               'W',    cell(2,1));
for i=1:numPoints
    [xData(1).II(i), xData(1).JJ(i)] = ind2sub([Nlon,Nlat], L2C_L_1D(i));
    [xData(2).II(i), xData(2).JJ(i)] = ind2sub([NX,NX*6], L2C_C_1D(i));
end
xData(1).W = L2C_S;
xData(2).W = L2C_S; % do not use!

%----------------------------------------------
% Correct xData to match GMAO tilefile mapping
%---------------------------------------------
% Correct longitude shift
%lonOffset = -2;
lonOffset = -10/(360/Nlon);
xData(1).II = shiftIndexes( xData(1).II, lonOffset, Nlon );

% Swap faces
faceRemapping =  [4 5 1 2 6 3];
xData(2) = swapCSFaces( xData(2), faceRemapping);

% Flip face 6 in both dimensions
xData(2) = flipCSFaces( xData(2), 6, '2d');

% Transpose and flip (II only) faces 3-5
xData(2) = transposeCSFaces( xData(2), [3 4 5]);
xData(2) = flipCSFaces( xData(2), [3 4 5], 'II');