function [dim] = pos2dim(pos)

% POS2DIM reconstructs the volumetric dimensions from an ordered list of 
% positions.
%
% Use as
%  [dim] = pos2dim(pos) where pos is an ordered list of positions
%  The output dim is a 3-element vector whichi correspond to the 3D 
%  volumetric dimensions

% Copyright (C) 2009, Jan-Mathijs Schoffelen

if isstruct(pos),
  %the input is a structure
  dimord = pos.dimord;
  dimtok = tokenize(dimord, '_');
  pos    = pos.pos;
end

%this part depends on the assumption that the list of positions is describing a full 3D volume in 
%an ordered way which allows for the extraction of a transformation matrix
%i.e. slice by slice
npos = size(pos,1);
dpos = zscore(abs(diff(pos,[],1)));

[tmp, ind] = max(dpos,[],2);
dim(1)     = find(tmp>1.5,1,'first');
dpos       = dpos(dim:dim:npos-1,:);
[tmp, ind] = max(dpos(:,setdiff(1:3, ind(dim))),[],2);
dim(2)     = find(tmp>1.1*min(tmp),1,'first'); %this threshold seems to work on what I tried out
dim(3)     = npos./prod(dim);
