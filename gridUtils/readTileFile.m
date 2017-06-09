function [ xData ] = readTileFile( tileFile )
%READTILEFILE Read a MAPL regrid tile file
%   Input is the path of the unformatted FORTRAN binary tile file. The
%   output, xData, is converted to a regridding object by the function
%   "genRegridObj".

assert(fastExist(tileFile),'readTileFile:fileNotFound','Tile file does not exist');
fID = fopen(tileFile,'rb');
cleanupFn = onCleanup(@()(fclose(fID)));

% Number of points
readFORTRANRecord(fID,'*int32',4);
nGrids = readFORTRANRecord(fID,'*int32',4);
assert(nGrids == 2,'readTileFile:badGridCount','Number of grids must be 2');

xData = struct();
for iGrid = 1:nGrids
    xData(iGrid).Name = strtrim(char(readFORTRANRecord(fID,'*char',1))');
    xData(iGrid).NX = readFORTRANRecord(fID,'*int32',4);
    xData(iGrid).NY = readFORTRANRecord(fID,'*int32',4);
end

% Type?
readFORTRANRecord(fID,'single',4);
% X locations?
readFORTRANRecord(fID,'single',4);
% Y locations?
readFORTRANRecord(fID,'single',4);

for iGrid = 1:nGrids
    % II and JJ
    xData(iGrid).II = round(readFORTRANRecord(fID,'single',4));
    xData(iGrid).JJ = round(readFORTRANRecord(fID,'single',4));

    % Weights
    xData(iGrid).W = readFORTRANRecord(fID,'single',4);
end
end

