function [ xData ] = displayTileFile( tileFile )
% EXPLORETILEFILE Read a MAPL regrid tile file for exploration/validation 
%    Input: Path of the unformatted FORTRAN binary tilefile
%    Output: xData struct array
%
% NOTES: Type, x/y locations, and a mystery data variable are not used by 
%        MAPL (as far as we can tell) but are included in xData as grid1 
%        values during this read for exploration and duplication purposes.
%
% Lizzie Lundgren, 9/26/16

assert(fastExist(tileFile),'readTileFile:fileNotFound', ...
       'Tile file does not exist');
fID = fopen(tileFile,'rb');
cleanupFn = onCleanup(@()(fclose(fID)));

fprintf('\n-----------------------------------------------\n');
fprintf('READING BINARY TILEFILE:\n');
fprintf('  %s\n', tileFile);

% Number of points
[pdata,rlen]=readFORTRANRecord(fID,'*int32',4);
fprintf('   # points: %d, mem = %d bytes\n', pdata, rlen);

% Number of grids
[nGrids, rlen] = readFORTRANRecord(fID,'*int32',4);
fprintf('   # grids: %d, mem = %d bytes\n', nGrids, rlen);
assert(nGrids == 2,'readTileFile:badGridCount', ...
       'Number of grids must be 2');

xData = struct();
for iGrid = 1:nGrids

    % Grid name	
    [name, rlen] = readFORTRANRecord(fID,'*char',1);	
    xData(iGrid).Name = strtrim(char(name)');
    trans_name = xData(iGrid).Name;
    fprintf('   Grid %d Name: %s, mem = %d bytes\n', ...
	    iGrid, trans_name, rlen);

    % Grid NX
    [nx,rlen] = readFORTRANRecord(fID,'*int32',4);
    xData(iGrid).NX = nx;
    fprintf('   Grid %d NX: %d, mem = %d bytes\n', ...
	    iGrid, nx, rlen);

    % Grid NY
    [ny,rlen] = readFORTRANRecord(fID,'*int32',4);
    xData(iGrid).NY = ny;
    fprintf('   Grid %d NY: %d, mem = %d bytes\n', ...
	    iGrid, ny, rlen);
end


% Type (not used)
[tdata, rlen] = readFORTRANRecord(fID,'single',4);
tlen = length(tdata);
fprintf('   Type: min = %d, max = %d, length = %d, mem = %d bytes\n', ...
       min(tdata), max(tdata), tlen, rlen);
xData(1).Type = tdata;

% X locations (not used)
[xdata, rlen] = readFORTRANRecord(fID,'single',4);
xlen = length(xdata);
fprintf('   Xloc: min = %d, max = %d, length = %d, mem = %d bytes\n', ...
       min(xdata), max(xdata), xlen, rlen);
xData(1).Xloc = xdata;

% Y locations (not used)
[ydata,rlen] = readFORTRANRecord(fID,'single',4);
ylen = length(ydata);
fprintf('   Yloc: min = %d, max = %d, length = %d, mem = %d bytes\n', ...
       min(ydata), max(ydata), ylen, rlen);
xData(1).Yloc = ydata;

for iGrid = 1:nGrids

    % Grid II
    [iidata, rlen] = readFORTRANRecord(fID,'single',4);
    meanii = mean(iidata(:));
    lenii = size(iidata,1);
    fprintf(['   Grid %d II: mean = %d, length = %d, ', ...
	    'mem = %d bytes\n'], iGrid, meanii, lenii, rlen);
    xData(iGrid).II = iidata;

    % Grid II
    [jjdata, rlen] = readFORTRANRecord(fID,'single',4);
    meanjj = mean(jjdata(:));
    lenjj = size(jjdata,1);
    fprintf(['   Grid %d JJ: mean = %d, length = %d, ', ...
	    'mem = %d bytes\n'], iGrid, meanjj, lenjj, rlen);
    xData(iGrid).JJ = jjdata;

    % Weights
    [wdata, rlen] = readFORTRANRecord(fID,'single',4);
    meanw = mean(wdata(:));
    lenw = size(wdata,1);
    fprintf(['   Grid %d W: mean = %d, length = %d, ', ...
	    'mem = %d bytes\n'], iGrid, meanw, lenw, rlen);
    xData(iGrid).W = wdata;
end

% Mystery variable (turn off by default)
%  NOTE: This variable exists in GMAO tilefiles but is not written
%        to new tilefiles with writeTileFile.m. For this reason, it
%        is commented out below. Uncomment this code to explore the
%        full GMAO tilefile data set.
%[testdata,rlen] = readFORTRANRecord(fID,'single',4);
%meantest = mean(testdata(:));
%lentest = size(testdata,1);
%xData(1).mysteryVar = testdata;
%fprintf(['   mysteryVar: mean = %d, length = %d, ', ...
%	'mem = %d bytes\n'], meantest, lentest, rlen);

end

