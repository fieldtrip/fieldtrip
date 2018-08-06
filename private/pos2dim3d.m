function [dim] = pos2dim3d(pos,dimold)

% POS2DIM3D reconstructs the volumetric dimensions from an ordered list of 
% positions. optionally, the original dim can be provided, and the (2:end)
% elements are appended to the output.
%
% Use as
%   [dim] = pos2dim3d(pos, dimold)
% where pos is an ordered list of positions and where the (optional)
% dimold is a vector with the original dimensionality of the anatomical
% or functional data.
%
% The output dim is a 1x3 or 1xN vector of which the first three elements
% correspond to the 3D volumetric dimensions.
%
% See also POS2DIM, POS2TRANSFORM

% Copyright (C) 2009, Jan-Mathijs Schoffelen

if nargin==1 && ~isstruct(pos),
  dimold = zeros(0,2);
elseif isstruct(pos),
  % the input is a FieldTrip data structure
  dimord = pos.dimord;
  dimtok = tokenize(dimord, '_');
  for i = 1:length(dimtok)
    if strcmp(dimtok{i},'pos'),
      dimold(i,1) = size(pos.pos,1);
    else
      dimold(i,1) = numel(getfield(pos, dimtok{i}));
    end
  end
  pos = pos.pos;
else
  if size(pos,1)~=dimold(1),
    ft_error('the first element in the second input should be equal to the number of positions');
  end
end

% extract the dim now that the bookkeeping is done
dim = pos2dim(pos);

