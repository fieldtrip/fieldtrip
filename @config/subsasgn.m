function y = subsasgn(x, index, val)

% SUBSASGN Assign a new value to a specified field in a config objects and increment its assignment counter.

% Copyright (C) 2012-2015, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if isa(val, 'config')
  warning('nested configs are not allowed, converting to struct');
  val = struct(val);
end

if length(index)==1
  switch index.type
    case '.'
      y = set(x, index.subs, val);
    case '{}'
      y = set(x, index.subs, val);
      % error('Cell contents reference from a non-cell-array object.');
    case '()'
      y = set(x, index.subs, val);
      % error('Index exceeds matrix dimensions.');
    otherwise
      error('Incorrect contents reference');
  end
elseif length(index)>1 && strcmp(index(1).type, '()')
  
  ndims = length(index(1).subs);
  
  %   if ndims==1 && isvector(x) && ~isequal(index(1).subs{1}, ':')
  %     index(1).subs = {1 index(1).subs{1}};
  %     ndims = 2;
  %   end
  
  % copy the array to the output variable
  y = x;
  
  if ndims==1
    % convert the array into a column vector
    y = reshape(y, numel(y), 1);
  end
  
  for i=1:ndims
    if isequal(index(1).subs{i}, ':')
      index(1).subs{i} = 1:size(x, i);
    end
  end
  
  % add singleton dimensions to the end
  index(1).subs((ndims+1):6) = {1};
  
  % the array may need to be extended in size
  for i=1:ndims
    dimsize = size(y, i);
    maxindx = max(index(1).subs{i});
    if dimsize<maxindx
      % construct an array like {'name1', [], 'name2', [], ...}
      fn = fieldnames(y);
      fv = repmat({[]}, size(fn));
      ff = cat(1, fn', fv');
      ff = config(ff{:});
      switch i
        case 1
          y((dimsize+1):maxindx,:,:,:,:,:) = ff;
        case 2
          y(:,(dimsize+1):maxindx,:,:,:,:) = ff;
        case 3
          y(:,:,(dimsize+1):maxindx,:,:,:) = ff;
        case 4
          y(:,:,:,(dimsize+1):maxindx,:,:) = ff;
        case 5
          y(:,:,:,:,(dimsize+1):maxindx,:) = ff;
        case 6
          y(:,:,:,:,:,(dimsize+1):maxindx) = ff;
        otherwise
          error('unsupported style of indexing');
      end % case
    end
  end
  
  if ~isfield(y, index(2).subs)
    % an empty field should be added to all elements of the array
    for i=1:numel(y)
      y(i) = subsasgn(y(i), index(2:end), []);
    end
  end
  
  % the value of the field should be changed for the specific elements of the array
  s1 = index(1).subs{1};
  s2 = index(1).subs{2};
  s3 = index(1).subs{3};
  s4 = index(1).subs{4};
  s5 = index(1).subs{5};
  s6 = index(1).subs{6};
  y(s1, s2, s3, s4, s5, s6) = subsasgn(y(s1, s2, s3, s4, s5, s6), index(2:end), val);
  
  if ndims==1
    if isscalar(x) && ~isscalar(y)
      y = reshape(y, 1, numel(y));
    elseif isrow(x) && ~isrow(y)
      y = reshape(y, 1, numel(y));
    elseif iscolumn(x) && ~iscolumn(y)
      y = reshape(y, numel(y), 1);
    elseif numel(x) == numel(y)
      y = reshape(y, size(x));
    elseif numel(y)>1
      error('unsupported style of indexing');
    end
  end
  
else
  % use recursion to find the subfield that is being indexed
  if ~isfield(x, index(1).subs)
    % the subfield does not exists, create a new config subfield
    switch index(2).type
      case '.'
        x  = subsasgn(x, index(1), config());
      case '()'
        x  = subsasgn(x, index(1), []);
      case '{}'
        x  = subsasgn(x, index(1), {});
      otherwise
        error('Incorrect contents reference');
    end
  end
  y1 = subsref(x, index(1), false);
  y2 = subsasgn(y1, index(2:end), val);
  y  = set(x, index(1).subs, y2);
end

