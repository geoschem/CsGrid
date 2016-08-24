function fEx = fastExist(fileName)
%FASTEXIST Check for existence of a file without searching MATLAB path

dirFile = dir(fileName);
if length(dirFile) == 1
   fEx = ~(dirFile.isdir);
else
   fEx = false;
end

end

