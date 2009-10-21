function [dim] = pos2dim3d(pos,dimold)

% POS2DIM3D reconstructs the volumetric dimensions from an ordered list of 
% positions. optionally, the original dim can be provided, and the (2:end)
% elements are appended to the output.
%
% Use as
%  [dim] = pos2dim3d(pos, dimold) where pos is an ordered list of positions
%  and dimold optionally the original dimensionality of the functional data
%  The output dim is a 3+ element vector of which the first three elements
%  correspond to the 3D volumetric dimensions

% Copyright (C) 2009, Jan-Mathijs Schoffelen
% $Log: pos2dim3d.m,v $
% Revision 1.5  2009/08/12 13:57:08  jansch
% added log for previous changes: changed standardise into zscore and allowed
% for single input argument as a source structure
%
% Revision 1.4  2009/08/12 13:56:06  jansch
% *** empty log message ***
%
% Revision 1.3  2009/07/01 10:53:34  jansch
% updated help and added log
%

if nargin==1 && ~isstruct(pos),
  dimold = zeros(0,2);
elseif isstruct(pos),
  %the input is a structure
  dimord = pos.dimord;
  dimtok = tokenize(dimord, '_');
  for i = 1:length(dimtok)
    if strcmp(dimtok{i},'pos'),
      dimold(i,1) = size(pos.pos,1);
    else
      dimold(i,1) = numel(getfield(pos, dimtok{i}));
    end
  end
  pos    = pos.pos;
else
  if size(pos,1)~=dimold(1),
    error('the first element in the second input should be equal to the number of positions');
  end
end

%this part depends on the assumption that the list of positions is describing a full 3D volume in 
%an ordered way which allows for the extraction of a transformation matrix
%i.e. slice by slice
npos = size(pos,1);
dpos = zscore(abs(diff(pos,[],1)));

[tmp, ind] = max(dpos,[],2);
tmpdim(1)  = find(tmp>2,1,'first');
dpos       = dpos(tmpdim:tmpdim:npos-1,:);
[tmp, ind] = max(dpos(:,setdiff(1:3, ind(tmpdim))),[],2);
tmpdim(2)  = find(tmp>1.1*min(tmp),1,'first'); %this threshold seems to work on what I tried out
tmpdim(3)  = npos./prod(tmpdim);
dim        = [tmpdim dimold(2:end)];
