function [ xDataOut ] = flipCSFaces( xDataIn, faceVec, faceDim )
% FLIPCSFACE Flip xData struct cubed sphere indexes along one dimension for
%            user-specified faces. 
%
%    Input:  (1) xDataIn:  xData struct (length 1), with xData.II and xData.JJ
%            (2) faceVec:  face number array with integers in set [1,6]
%            (3) faceDim:  string indicating dimension to flip, 'II', 'JJ', 
%                          or '2d' (flips both II and JJ)
%
%    Output: (1) xDataOut: xData struct (length 1), with faceDim indexes flipped
%
%    Use Case: Flip xData(2).II (cs) indexes on the 6th face: 
%              >> xData(2).II = flipCSFaces(xData(2).II, 6, 'II');
%
% Lizzie Lundgren, 9/29/16

% Initialize
xDataOut = xDataIn;
f = 1;

% Check that inputs are valid (in development)
validDims={'II' 'JJ' '2d'};
if ~any(ismember(validDims, faceDim))
  error('invalid faceDim passed to flipCSFaces.m');
end

% Loop over faces
while f <= length(faceVec);
  faceNum = faceVec(f);
  f = f+1;

  % Get cubed sphere side length
  NX = xDataIn.NX; 

  % Set min/max indexes based on face number and cubed sphere side length
  minJJind = NX*(faceNum-1)+1;
  maxJJind = NX*faceNum;
  minIIind = 1;
  maxIIind = NX;
  
  % Find elements for this face based on xData.JJ values (NX*6 side)
  thisFaceInd = find( xDataIn.JJ >= minJJind & xDataIn.JJ <= maxJJind );
  
  % Loop over the elements found, flipping the user-specified dimension
  for i = 1:length(thisFaceInd);
    if strcmp(faceDim, 'JJ') | strcmp(faceDim, '2d');
          xDataOut.JJ(thisFaceInd(i)) = minJJind + ( maxJJind ...
				   - xDataIn.JJ(thisFaceInd(i)) );
    end
    if strcmp(faceDim, 'II') | strcmp(faceDim, '2d');
          xDataOut.II(thisFaceInd(i)) = minIIind + ( maxIIind ...
				   - xDataIn.II(thisFaceInd(i)) );
    end
  end 
end
