function [ USort, VSort ] = resortUV( U, V )
%RESORTUV Re-orients U and V when read from GCHP data

assert(numel(U) == numel(V) && ndims(U) == ndims(V) && ...
    all(size(U) == size(V)),'resortUV:inconsistentInput',...
    'Inputs U and V must be identically sized')

if ismatrix(U)
    warning('resortUV:reorderingData','Reordering data from 2D to 3D - assuming GCHP-sourced.');
    U = reorderCS(U,true);
    V = reorderCS(V,true);
end

USort = U;
VSort = V;

% Faces
% 1: "Side" - Greenwich
% 2: "Side" - East
% 3: "South pole"
% 4: "Side" - Crossing date line
% 5: "Side" - West
% 6: "North pole"
%
%          6
% -------------------
%  4 | 5 | 1 | 2 | 4 
% -------------------
%          3
%
USrc = ([ 1, 0;...
          1, 0;...
          0,-1;...
         -1, 0;...
         -1, 0;...
          0,-1]);
VSrc = ([ 0, 1;...
          0, 1;...
          1, 0;...
          0,-1;...
          0,-1;...
          1, 0]);
      
for iFace = 1:6
    USort(:,:,iFace) = U(:,:,iFace).*USrc(iFace,1) + V(:,:,iFace).*USrc(iFace,2);
    VSort(:,:,iFace) = U(:,:,iFace).*VSrc(iFace,1) + V(:,:,iFace).*VSrc(iFace,2);
end

end

