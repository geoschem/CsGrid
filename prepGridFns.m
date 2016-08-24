function [  ] = prepGridFns(  )
%PREPGRIDFNS Prepare the path

dirList = {'GridData','gridUtils'};
for iDir = 1:length(dirList)
    currDir = dirList{iDir};
    assert(exist(currDir,'dir')==7,'prepGridFns:missingFolder','Folder %s not present',currDir);
end
addpath(genpath('gridUtils'));


end

