function [ indArrOut ] = shiftIndexes( indArrIn, indOffset, indMax )
% SHIFTINDEXES Apply an offset to an array of indexes.
%
%    Inputs:  (1) indArrIn:  array of input indexes (integers)
%             (2) offset:    positive or negative integer offset
%             (3) indMax:    maximum allowed index in the array
%
%    Outputs: (1) indArrOut: array of input indexes shifted by offset.
%                            Indexes wrap-around so as not to go below 
%                            1 or above indMax following offset application.
%    
%    NOTE: While there is user-specified maximum index, the minimum is
%          assumed to be 1.
 
% Lizzie Lundgren

indArrOut = zeros(length(indArrIn),1);
for i = 1:length(indArrIn);
    if indOffset > 0 &  indArrIn(i) > (Nlon-indOffset);
        indArrOut(i) = indArrIn(i) + indOffset - Nlon;
    elseif indOffset < 0 & indArrIn(i) < abs(indOffset)+1; 
        indArrOut(i) = indArrIn(i) + indOffset + Nlon;
    else
        indArrOut(i) = indArrIn(i) + indOffset;
    end
end
