function [ xDataOut ] = transposeCSFaces( xDataIn, faceVec )
% SWAPCSFACES Transpose xData struct cubed sphere indexes for each face in
%             a vector of user-specified faces.
%
%    Input:  (1) xDataIn:  xData struct (length 1), with xData.II and xData.JJ
%            (2) faceVec:  face number array with integers in set [1,6]
%
%    Output: (1) xDataOut: xData struct (length 1), with indexes II and JJ
%                          of faces in faceVec array transposed
%
%    Use Case: Transpose xData(2) indexes for faces 3, 4, and 5.
%              >> xData(2) = transposeCSFaces( xData(2), [3 4 5]);
%
% Lizzie Lundgren, 9/29/16

% Initialize
xDataOut = xDataIn;
NX = xDataIn.NX; 
f = 1;

% Check that inputs are valid
assert( min(faceVec) >= 1 | max(faceVec) <= 6, ...
       'invalid face number in faceVec passed to transposeCSFaces.m');
assert( isfield(xDataIn,'II'), ...
       'II field does not exist in struct passed to transposeCSFaces.m');
assert( isfield(xDataIn,'JJ'), ...
       'JJ field does not exist in struct passed to transposeCSFaces.m');
assert( isfield(xDataIn,'NX'), ...
       'NX field does not exist in struct passed to transposeCSFaces.m');

% Loop over faces
while f <= length(faceVec);
  faceNum = faceVec(f);
  f = f+1;

  % Set min/max indexes based on face number and cubed sphere side length
  minJJind = NX*(faceNum-1)+1;
  maxJJind = NX*faceNum;
  minIIind = 1;
  maxIIind = NX;
  
  % Find elements for this face based on xData.JJ values (NX*6 side)
  thisFaceInd = find( xDataIn.JJ >= minJJind & xDataIn.JJ <= maxJJind );
  
  % Loop over the elements found
  for i = 1:length(thisFaceInd);
     xDataOut.II(thisFaceInd(i)) = xDataIn.JJ(thisFaceInd(i)) - minJJind + 1;
     xDataOut.JJ(thisFaceInd(i)) = xDataIn.II(thisFaceInd(i)) + minJJind - 1;
  end 
end
