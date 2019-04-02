function varargout = subsref(x, index, inc)

% SUBSREF Return the value of a specified field in a config objects and increment its reference counter.

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

if nargin<3
  inc = true;
end

if strcmp(index(1).type, '()')
  ndims = numel(index(1).subs);
  if ndims==1
    % convert the array into a column vector
    x = reshape(x, numel(x), 1);
  end
  for i=1:ndims
    if isequal(index(1).subs{i}, ':')
      index(1).subs{i} = 1:size(x, i);
    end
  end
end

if length(index)==1
  switch index.type
    case '.'
      varargout      = cell(1, numel(x));
      [varargout{:}] = get(x, index.subs, inc);
    case '{}'
      error('Cell contents reference from a non-cell-array object.');
    case '()'
      % add singleton dimensions to the end
      index.subs((ndims+1):6) = {1};
      s1 = index.subs{1};
      s2 = index.subs{2};
      s3 = index.subs{3};
      s4 = index.subs{4};
      s5 = index.subs{5};
      s6 = index.subs{6};
      % get the selection
      varargout{1} = x(s1, s2, s3, s4, s5, s6);
    otherwise
      error('Incorrect contents reference');
  end
else
  % use recursion to find the subfield that is being indexed
  if ischar(index(1).subs)
    % this happens when indexing a single structure
    varargout = cell(1, 1);
  elseif iscell(index(1).subs)
    % this happens when indexing an array
    varargout = cell(1, prod(cellfun(@numel, index(1).subs)));
  end
  try
    [varargout{:}] = subsref(subsref(x, index(1)), index(2:end));
  catch
    me = lasterror;
    if strcmp(me.identifier,  'MATLAB:unassignedOutputs')
      error('Scalar index required for this type of multi-level indexing.');
    else
      error(me.message);
    end
  end
end

% TEST CODE
% function varargout = subsref(x, index, inc)
% % SUBSREF Return the value of a specified field in a config objects and increment its reference counter.
%
% y = [];
% if nargin<3
%   inc = true;
% end
%
% if length(index)==1
%   switch index.type
%     case '.'
%       if numel(x)>1
%         y = cell(1,numel(x));
%         % fields from multiple (sub)structures are requested, loop over each
%         for iobj = 1:numel(x)
%           y{iobj} = get(x(iobj), index.subs, inc);
%         end
%       else
%         y{1} = get(x, index.subs, inc);
%       end
%     case '{}'
%       error('Cell contents reference from a non-cell-array object.');
%     case '()'
%         y{1} = x(index.subs{1});
%     otherwise
%       error('Incorrect contents reference');
%   end
% else
%   % use recursion to find the subfield that is being indexed
%     y{1} = subsref(subsref(x, index(1)), index(2:end));
% end
% varargout = y;
