function [ xData ] = exploreTileFile( tileFile )
% EXPLORETILEFILE Read a MAPL regrid tile file for exploration/validation 
%    Input: Path of the unformatted FORTRAN binary tilefile
%    Output: xData struct array
%
% NOTES: Type, x/y locations, and a mystery data variable are not used by 
%        MAPL (as far as we can tell) but are included in xData as grid1 
%        values during this read for exploration and duplication purposes.
%
% Lizzie Lundgren, 9/26/16

assert(fastExist(tileFile),'readTileFile:fileNotFound','Tile file does not exist');
fID = fopen(tileFile,'rb');
cleanupFn = onCleanup(@()(fclose(fID)));

fprintf('\n-----------------------------------------------\n');
fprintf('READING BINARY TILEFILE:\n');
fprintf('  %s\n', tileFile);

% Number of points
[pdata,rlen]=readFORTRANRecord(fID,'*int32',4);
fprintf('   # points = %d, mem = %d bytes\n', pdata, rlen);

% Number of grids
[nGrids, rlen] = readFORTRANRecord(fID,'*int32',4);
fprintf('   # grids = %d, mem = %d bytes\n', nGrids, rlen);
assert(nGrids == 2,'readTileFile:badGridCount','Number of grids must be 2');

xData = struct();
for iGrid = 1:nGrids

    % Grid name	
    [name, rlen] = readFORTRANRecord(fID,'*char',1);	
    xData(iGrid).Name = strtrim(char(name)');
    trans_name = xData(iGrid).Name;
    fprintf('   Grid %d Name = %s, mem = %d bytes\n', iGrid, trans_name, rlen);

    % Grid NX
    [nx,rlen] = readFORTRANRecord(fID,'*int32',4);
    xData(iGrid).NX = nx;
    fprintf('   Grid %d NX = %d, mem = %d bytes\n', iGrid, nx, rlen);

    % Grid NY
    [ny,rlen] = readFORTRANRecord(fID,'*int32',4);
    xData(iGrid).NY = ny;
    fprintf('   Grid %d NY = %d, mem = %d bytes\n', iGrid, ny, rlen);
end


% Type
[tdata, rlen] = readFORTRANRecord(fID,'single',4);
xData(1).Type = tdata;

% X locations
[xdata, rlen] = readFORTRANRecord(fID,'single',4);
xData(1).Xloc = xdata;

% Y locations
[ydata,rlen] = readFORTRANRecord(fID,'single',4);
xData(1).Yloc = ydata;

for iGrid = 1:nGrids

    % Grid II
    [iidata, rlen] = readFORTRANRecord(fID,'single',4);
    meanii = mean(iidata(:));
    lenii = size(iidata,1);
    fprintf('   Grid %d mean II = %d, length II = %d, mem = %d bytes\n', iGrid, meanii, lenii, rlen);
    xData(iGrid).II = round(iidata);

    % Grid II
    [jjdata, rlen] = readFORTRANRecord(fID,'single',4);
    meanjj = mean(jjdata(:));
    lenjj = size(jjdata,1);
    fprintf('   Grid %d mean JJ = %d, length JJ = %d, mem = %d bytes\n', iGrid, meanjj, lenjj, rlen);
    xData(iGrid).JJ = round(jjdata);

    % Weights
    [wdata, rlen] = readFORTRANRecord(fID,'single',4);
    meanw = mean(wdata(:));
    lenw = size(wdata,1);
    fprintf('   Grid %d mean W = %d, length W = %d, mem = %d bytes\n', iGrid, meanw, lenw, rlen);
    xData(iGrid).W = wdata;
end

% Mystery variable
[testdata,rlen] = readFORTRANRecord(fID,'single',4);
meantest = mean(testdata(:));
lentest = size(testdata,1);
xData(1).mysteryVar = testdata;
fprintf('   mysteryVar mean = %d, length mysteryVar = %d, mem = %d bytes\n', meantest, lentest, rlen);

end

