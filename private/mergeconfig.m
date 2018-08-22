function original = mergeconfig(original, default, emptymeaningful)

% MERGECONFIG merges the fields of a structure with defaults into an original
% configuration structure. The default is only copied in case the field is absent in
% the original.
%
% Use as
%   output = mergeconfig(original, default)

% Copyright (C) 2009-2012, Robert Oostenveld
%
% $Id$

if nargin<3
  emptymeaningful = true;
end

if isempty(original) && ~isstruct(original)
  % ensure that it is an empty struct, not an empty double
  original = struct;
end

if isempty(default) && ~isstruct(default)
  % ensure that it is an empty struct, not an empty double
  default = struct;
end

if ~emptymeaningful
  original = remove_empty(original);
  default = remove_empty(default);
end

% merge the original with the fields from default
defaultfields = fieldnames(default);
inputfields = fieldnames(original);
allfields = union(defaultfields, inputfields);

if numel(default)>1 && numel(default)==numel(original)
  % create an empty structure that has all the fields
  tmp = empty_struct(allfields);
  for i=1:numel(original)
    for j=1:numel(allfields)
      if isfield(original(i), allfields{j})
        tmp(i).(allfields{j}) = original(i).(allfields{j});
      else
        tmp(i).(allfields{j}) = default(i).(allfields{j});
      end
    end
  end
  original = tmp;
  
elseif numel(default)~=numel(original)
  for j=1:numel(allfields)
    % ensure that the original has all fields
    if ~isfield(original, allfields{j})
      original(1).(allfields{j}) = [];
    end
    % ensure that the default has all fields
    if ~isfield(default, allfields{j})
      default(1).(allfields{j}) = [];
    end
  end
  % simply concatenate them
  original = cat(1, original(:), default(:));
  
else
  for i=1:length(defaultfields)
    fn = defaultfields{i};
    if  ~isfield(original, fn) && ~isstruct(default.(fn))
      % simply copy the value over
      original.(fn) = default.(fn);
    elseif ~isfield(original, fn) &&  isstruct(default.(fn))
      % simply copy the substructure over
      original.(fn) = default.(fn);
    elseif  isfield(original, fn) && ~isstruct(default.(fn))
      % do not copy it over, keep the original value
    elseif  isfield(original, fn) &&  isstruct(default.(fn))
      % merge the two substructures using recursive call
      original.(fn) = mergeconfig(original.(fn), default.(fn));
    end
  end % for all default fields
  
end % dealing with struct-arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = remove_empty(s)
if isempty(s)
  return
else
  fn = fieldnames(s);
  fn = fn(structfun(@isempty, s));
  s = removefields(s, fn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = empty_struct(fn)
s = struct;
for i=1:numel(fn)
  s.(fn{i}) = [];
end
s = s([]);