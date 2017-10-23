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
defaultfields = fieldnames(default);
inputfields = fieldnames(input);
allfields = union(defaultfields, inputfields);

if numel(default)>1 && numel(default)==numel(input)
  % create an empty structure that has all the fields
  tmp = emptystruct(allfields);
  for i=1:numel(input)
    for j=1:numel(allfields)
      if isfield(input(i), allfields{j})
        tmp(i).(allfields{j}) = input(i).(allfields{j});
      else
        tmp(i).(allfields{j}) = default(i).(allfields{j});
      end
    end
  end
  input = tmp;
  
elseif numel(default)~=numel(input)
  for j=1:numel(allfields)
    % ensure that the input has all fields
    if ~isfield(input, allfields{j})
      input(1).(allfields{j}) = [];
    end
    % ensure that the default has all fields
    if ~isfield(default, allfields{j})
      default(1).(allfields{j}) = [];
    end
  end
  % simply concatenate them
  input = cat(1, input(:), default(:));
  
else
  for i=1:length(defaultfields)
    fn = defaultfields{i};
    if  ~isfield(input, fn) && ~isstruct(default.(fn))
      % simply copy the value over
      input.(fn) = default.(fn);
    elseif ~isfield(input, fn) &&  isstruct(default.(fn))
      % simply copy the substructure over
      input.(fn) = default.(fn);
    elseif  isfield(input, fn) && ~isstruct(default.(fn))
      % do not copy it over, keep the original value
    elseif  isfield(input, fn) &&  isstruct(default.(fn))
      % merge the two substructures using recursive call
      input.(fn) = mergeconfig(input.(fn), default.(fn));
    end
  end % for all default fields
  
end % dealing with struct-arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = emptystruct(fn)
s = struct;
for i=1:numel(fn)
  s.(fn{i}) = [];
end
s = s([]);