function s1 = mergestruct(s1, s2, emptymeaningful)

% MERGESTRUCT merges the fields of a structure with another structure. The fields of
% the 2nd structure are only copied in case they are absent in the 1st structure.
%
% Use as
%   s3 = mergestruct(s1, s2, emptymeaningful)
%
% See also PRINTSTRUCT, APPENDSTRUCT, COPYFIELDS, KEEPFIELDS, REMOVEFIELDS, MERGETABLE

% Copyright (C) 2009-2022, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<3
  emptymeaningful = true;
end

if isempty(s1) && ~isstruct(s1)
  % ensure that it is an empty struct, not an empty double
  s1 = struct;
end

if isempty(s2) && ~isstruct(s2)
  % ensure that it is an empty struct, not an empty double
  s2 = struct;
end

if ~emptymeaningful
  s1 = remove_empty(s1);
  s2 = remove_empty(s2);
end

% merge the original with the fields from default
defaultfields = fieldnames(s2);
inputfields = fieldnames(s1);
allfields = union(defaultfields, inputfields);

if numel(s2)>1 && numel(s2)==numel(s1)
  % create an empty structure that has all the fields
  tmp = empty_struct(allfields);
  for i=1:numel(s1)
    for j=1:numel(allfields)
      if isfield(s1(i), allfields{j})
        tmp(i).(allfields{j}) = s1(i).(allfields{j});
      else
        tmp(i).(allfields{j}) = s2(i).(allfields{j});
      end
    end
  end
  s1 = tmp;
  
elseif numel(s2)~=numel(s1)
  for j=1:numel(allfields)
    % ensure that the original has all fields
    if ~isfield(s1, allfields{j})
      s1(1).(allfields{j}) = [];
    end
    % ensure that the default has all fields
    if ~isfield(s2, allfields{j})
      s2(1).(allfields{j}) = [];
    end
  end
  % simply concatenate them
  s1 = cat(1, s1(:), s2(:));
  
else
  for i=1:length(defaultfields)
    fn = defaultfields{i};
    if  ~isfield(s1, fn) && ~isstruct(s2.(fn))
      % simply copy the value over
      s1.(fn) = s2.(fn);
    elseif ~isfield(s1, fn) &&  isstruct(s2.(fn))
      % simply copy the substructure over
      s1.(fn) = s2.(fn);
    elseif  isfield(s1, fn) && ~isstruct(s2.(fn))
      % do not copy it over, keep the original value
    elseif  isfield(s1, fn) &&  isstruct(s2.(fn))
      % merge the two substructures using recursive call
      s1.(fn) = mergestruct(s1.(fn), s2.(fn));
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