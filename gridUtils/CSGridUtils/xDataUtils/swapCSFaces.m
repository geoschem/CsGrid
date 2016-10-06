function [ xDataOut ] = swapCSFaces( xDataIn, faceVec )
% SWAPCSFACES Swap cubed sphere faces in xData struct
%
%    Input:  (1) xDataIn:  xData struct (length 1), with xData.JJ field, the
%                          mapping indexes in set 1:NX*6
%            (2) faceVec:  array of new face indexes (length 6); array indexes
%                          are input faces, and array values are output faces
%
%    Output: (1) xDataOut: xData struct (length 1), with indexes JJ of faces
%                          in faceVec switched from face = faceVec index to
%                          face = faceVec value.
%
%    Use Case: Change xData(2).JJ indexes to swap faces: 1->4, 2->5, 3->1,
%              4->2, 5->6, and 6->3. This is the face mapping to go from 
%              Tempest to GMAO for tilefiles generation.
%              >> xData(2).JJ = ...
%                       swapCSFaces( xData(2).JJ, [4 5 1 2 6 3], NX);
%
% Lizzie Lundgren, 9/29/16

% Initialize
xDataOut = xDataIn;
NX = xDataIn.NX; 
numFaces = 6;

% Check that inputs are valid
assert( length(faceVec) == numFaces & length(unique(faceVec)) == numFaces, ...
        'invalid face mapping passed to swapCSFaces.m');
assert( isfield(xDataIn,'JJ'), ...
       'JJ field does not exist in struct passed to transposeCSFaces.m');
assert( isfield(xDataIn,'NX'), ...
       'NX field does not exist in struct passed to transposeCSFaces.m');

% Loop over faces
for f=1:numFaces;

   % define out blocks for faces (f is original face, a is new face)
   a=faceVec(f);

   % define the face block start and end indexes
   minJJind_in = NX*(f-1)+1;
   minJJind_out = NX*(a-1)+1;
   maxJJind_in = NX*f;

   % get JJ indexes corresponding to this face
   thisFaceInd = find(xDataIn.JJ >= minJJind_in & ...
	              xDataIn.JJ <= maxJJind_in);

   % Loop over all JJ indexes for this face		       
   for j = 1:length(thisFaceInd);

       % determine the offset of this JJ value from original face block start
       offsetJJ = xDataIn.JJ(thisFaceInd(j)) - minJJind_in;

       % set a new JJ value that is the same offset from new face block start
       xDataOut.JJ(thisFaceInd(j)) = minJJind_out + offsetJJ;
   end 
end

% Check that output is valid
assert(min(xDataOut.JJ) > 0, 'Invalid index array in swapFaces.');