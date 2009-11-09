function [type] = voltype(vol, desired)

% VOLTYPE determines the type of volume conduction model
%
% Use as
%   [type] = voltype(vol)
% to get a string describing the type, or
%   [flag] = voltype(vol, desired)
% to get a boolean value.
%
% See also READ_VOL, COMPUTE_LEADFIELD

% Copyright (C) 2007-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

if iscell(vol) && numel(vol)<4
  % this might represent combined EEG, ECoG and/or MEG
  type = cell(size(vol));
  if nargin<2
    desired = cell(size(vol)); % empty elements
  end
  for i=1:numel(vol)
    type{i} = voltype(vol{i}, desired{i});
  end
  return
end

if nargin<2
  % ensure that all input arguments are defined
  desired = [];
end

current_argin = {vol, desired};
if isequal(current_argin, previous_argin)
  % don't do the type detection again, but return the previous values from cache
  type = previous_argout{1};
  return
end

if isfield(vol, 'type')
  % preferably the structure specifies its own type
  type = vol.type;

elseif isfield(vol, 'r') && numel(vol.r)==1 && ~isfield(vol, 'label')
  type = 'singlesphere';

elseif isfield(vol, 'r') && isfield(vol, 'o') && isfield(vol, 'label')
  % this is before the spheres have been assigned to the coils
  % and every sphere is still associated with a channel
  type = 'multisphere';

elseif isfield(vol, 'r') && isfield(vol, 'o') && size(vol.r,1)==size(vol.o,1) && size(vol.r,1)>4
  % this is after the spheres have been assigned to the coils
  % note that this one is easy to confuse with the concentric one
  type = 'multisphere';

elseif isfield(vol, 'r') && numel(vol.r)>=2 && ~isfield(vol, 'label')
  type = 'concentric';

elseif isfield(vol, 'bnd')
  type = 'bem';

elseif isempty(vol)
  type = 'infinite';

else
  type = 'unknown';

end % if isfield(vol, 'type')

if ~isempty(desired)
  % return a boolean flag
  switch desired
    case 'bem'
      type = any(strcmp(type, {'bem', 'dipoli', 'asa', 'avo', 'bemcp', 'openmeeg'}));
    otherwise
      type = any(strcmp(type, desired));
  end
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {type};
if isempty(previous_argin)
  previous_argin  = current_argin;
  previous_argout = current_argout;
end

return % voltype main()
