function [ indArrOut ] = swapCSFaces( indArrIn, faceVec, NX )
% SWAPCSFACES Switch cubed sphere faces by switching the NX*NX arrays
%             within the cubed sphere array.
%
%    Input:  (1) indArrIn: array of cubed sphere indexes (dim w/ length NX*6)
%            (2) faceVec:  array of new face indexes (length 6); array value
%                          is new face which array index is original face
%            (3) NX:       length of cubed sphere face
%
%    Output: (1) indArrOut: array of cubed sphere indexes (dim w/ length NX*6)
%
% NOTES: Faces may need to be transposed or rotated following the swap.
%        Currently this implementation does not predict when that will need
%        to happen and does not automatically do it for you. This might
%        be possible to do, however.
%
% Lizzie Lundgren

% initialize output
indArrOut = zeros(length(indArrOut),1);

% loop over cubed sphere faces
for f=1:6;

   % define new blocks for faces (f is original face, a is new face)
   a=faceVec(f);

   % define the face block start and end indexes
   face_start_ind_orig = NX*(f-1)+1;
   face_end_ind_orig = NX*f;
   face_start_ind_new = NX*(a-1)+1;

   % get JJ indexes corresponding to this face
   this_face_ind = find(indArrIn >= face_start_ind_orig & ...
			indArrIn <= face_end_ind_orig);

   % Loop over all JJ indexes for this face		       
   for j = 1:length(this_face_ind);

       % determine the offset of this JJ value from original face block start
       offset_face_JJ = indArrIn(this_face_ind(j)) - face_start_ind_orig;

       % set a new JJ value that is the same offset from new face block start
       indArrOut(this_face_ind(j)) = face_start_ind_new + offset_face_JJ;
   end 
end

% Check that output is valid
assert(min(indArrOut) > 0, "Invalid index array in swapFaces.");