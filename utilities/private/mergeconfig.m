function input = mergeconfig(input, default)

% MERGECONFIG merges the fields of a default structure into an input
% configuration structure
%
% Use as
%   output = mergeconfig(input, default)

% Copyright (C) 2009-2012, Robert Oostenveld
%
% $Id$
  input = mergeconfig_helper(input, default);

function input=mergeconfig_helper(input,default)

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
  
  if isstruct(default) && numel(default)>1
    if isempty(input) || isstruct(input) && isempty(fieldnames(input))
      input = default;
    else
      for j=1:numel(default)
        if numel(input)<j
          input(j) = default(j);
        else
          input(j) = mergeconfig_helper(input(j), default(j));
        end
      end
    end
  elseif  ~isfield(input, fn) && ~isstruct(default.(fn))
    % simply copy the value over
    input.(fn) = default.(fn);
  elseif ~isfield(input, fn) &&  isstruct(default.(fn))
    % simply copy the substructure over
    input.(fn) = default.(fn);
  elseif  isfield(input, fn) && ~isstruct(default.(fn))
    % do not copy it over, keep the original value
  elseif  isfield(input, fn) &&  isstruct(default.(fn))
    % merge the two substructures using recursive call
    input.(fn) = mergeconfig_helper(input.(fn), default.(fn));
  end
  
end % for all fields
