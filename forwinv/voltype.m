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
% $Log: voltype.m,v $
% Revision 1.4  2009/04/01 06:51:43  roboos
% implemented caching for the type detection in case the same input is given multiple times
% this uses a persistent variable
%
% Revision 1.3  2009/02/02 13:05:22  roboos
% small fix in cpbem->bemcp
%
% Revision 1.2  2009/02/02 12:58:27  roboos
% added bemcp implementation
%
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.7  2008/04/30 13:40:34  roboos
% improved detection concentric eeg
%
% Revision 1.6  2008/04/16 08:04:33  roboos
% be flexible in determining whether it is bem
%
% Revision 1.5  2008/04/14 19:31:05  roboos
% return 'unknown' if the type cannot be determined
%
% Revision 1.4  2008/04/10 10:59:38  roboos
% better detection of multisphere meg (after preparing)
%
% Revision 1.3  2008/03/18 12:39:24  roboos
% change in comment, nothing functional
%
% Revision 1.2  2008/03/05 15:24:44  roboos
% changed the detection of various spherical models, added infinite vacuum
%
% Revision 1.1  2007/07/25 08:31:12  roboos
% implemented new helper function
%

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

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
      type = any(strcmp(type, {'bem', 'dipoli', 'asa', 'avo', 'bemcp'}));
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
