function input = mergeconfig(input, default)

% MERGECONFIG merges the fields of a default structure into an input
% configuration structure
%
% Use as
%   output = mergeconfig(input, default)

% Copyright (C) 2009-2012, Robert Oostenveld
%
% $Id$

if isempty(input) && ~isstruct(input)
  % ensure that it is an empty struct, not an empty double
  input = struct;
end

if isempty(default) && ~isstruct(default)
  % ensure that it is an empty struct, not an empty double
  default = struct;
end

% merge the input with the fields from default
fieldsused = fieldnames(default);
for i=1:length(fieldsused)
  fn = fieldsused{i};
  
  if     ~isfield(input, fn) && ~isstruct(default.(fn))
    % simply copy the value over
    input.(fn) = default.(fn);
  elseif ~isfield(input, fn) &&  isstruct(default.(fn))
    % simply copy the substructure over
    input.(fn) = default.(fn);
  elseif  isfield(input, fn) && ~isstruct(default.(fn))
    % do not copy it over, keep the original value
  elseif  isfield(input, fn) &&  isstruct(default.(fn))
    % merge the two substructures
    input.(fn) = mergeconfig(input.(fn), default.(fn));
  end
  
end % for all fields
