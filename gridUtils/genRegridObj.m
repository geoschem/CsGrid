function [xRegrid] = genRegridObj(xData)
%GENREGRIDOBJ Converts Tempest tile file data to a regridding array
% This is really just undoing some of the work done in
% readTempestGeneric, but it results in a standard interface so it's
% worthwhile (for now)
nEl = cast([xData(1).NX*xData(1).NY,xData(2).NX*xData(2).NY],'double');

% Much faster to call on a vector
fromVec = sub2ind([xData(1).NX,xData(1).NY],xData(1).II,xData(1).JJ);
toVec = sub2ind([xData(2).NX,xData(2).NY],xData(2).II,xData(2).JJ);

% Condense the information into a single sparse array
xRegrid = sparse(fromVec,toVec,xData(1).W,nEl(1),nEl(2));
end
