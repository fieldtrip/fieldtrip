function opt = ft_checkopt(opt, key, allowedtype, allowedval)

% FT_CHECKOPT does a validity test on the types and values of a configuration
% structure or cell-array with key-value pairs.
%
% Use as
%   opt = ft_checkopt(opt, key)
%   opt = ft_checkopt(opt, key, allowedtype)
%   opt = ft_checkopt(opt, key, allowedtype, allowedval)
%
% For allowedtype you can specify a string or a cell-array with multiple
% strings. All the default MATLAB types can be specified, such as
%   'double'
%   'logical'
%   'char'
%   'single'
%   'float'
%   'int16'
%   'cell'
%   'struct'
%   'function_handle'
% Furthermore, the following custom types can be specified
%   'doublescalar'
%   'doublevector'
%   'doublebivector'             i.e. [1 1] or [1 2]
%   'ascendingdoublevector'      i.e. [1 2 3 4 5], but not [1 3 2 4 5]
%   'ascendingdoublebivector'    i.e. [1 2], but not [2 1]
%   'doublematrix'
%   'numericscalar'
%   'numericvector'
%   'numericmatrix'
%   'charcell'
%
% For allowedval you can specify a single value or a cell-array
% with multiple values.
%
% This function will give an error or it returns the input configuration
% structure or cell-array without modifications. A match on any of the
% allowed types and any of the allowed values is sufficient to let this
% function pass.
%
% See also FT_GETOPT, FT_SETOPT

% Copyright (C) 2011-2012, Robert Oostenveld
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
  allowedtype = {};
end

if ~iscell(allowedtype)
  allowedtype = {allowedtype};
end

if nargin<4
  allowedval = {};
end

if ~iscell(allowedval)
  allowedval = {allowedval};
end

% get the value that belongs to this key
val  = ft_getopt(opt, key);       % the default will be []

if isempty(val) && ~any(strcmp(allowedtype, 'empty'))
  if isnan(ft_getopt(opt, key, nan))
    error('the option "%s" was not specified or was empty', key);
  end
end

% check that the type of the option is allowed
ok = isempty(allowedtype);
for i=1:length(allowedtype)
  switch allowedtype{i}
    case 'empty'
      ok = isempty(val);
    case 'charcell'
      ok = isa(val, 'cell') && all(cellfun(@ischar, val(:)));
    case 'doublescalar'
      ok = isa(val, 'double') && numel(val)==1;
    case 'doublevector'
      ok = isa(val, 'double') && sum(size(val)>1)==1;
    case 'ascendingdoublevector'
      ok = isa(val,'double') && issorted(val);
    case 'doublebivector'
      ok = isa(val,'double') && sum(size(val)>1)==1 && length(val)==2;
    case 'ascendingdoublebivector'
      ok = isa(val,'double') && sum(size(val)>1)==1 && length(val)==2 && val(2)>val(1);
    case 'doublematrix'
      ok = isa(val, 'double') && sum(size(val)>1)>1;
    case 'numericscalar'
      ok = isnumeric(val) && numel(val)==1;
    case 'numericvector'
      ok = isnumeric(val) && sum(size(val)>1)==1;
    case 'numericmatrix'
      ok = isnumeric(val) && sum(size(val)>1)>1;
    otherwise
      ok = isa(val, allowedtype{i});
  end
  if ok
    % no reason to do additional checks
    break
  end
end % for allowedtype

% construct a string that describes the type of the input variable
if isnumeric(val) && numel(val)==1
  valtype = sprintf('%s scalar', class(val));
elseif isnumeric(val) && numel(val)==length(val)
  valtype = sprintf('%s vector', class(val));
elseif isnumeric(val) && length(size(val))==2
  valtype = sprintf('%s matrix', class(val));
elseif isnumeric(val)
  valtype = sprintf('%s array', class(val));
else
  valtype = class(val);
end

if ~ok
  if length(allowedtype)==1
    error('the type of the option "%s" is invalid, it should be "%s" instead of "%s"', key, allowedtype{1}, valtype);
  else
    error('the type of the option "%s" is invalid, it should be any of %s instead of "%s"', key, printcell(allowedtype), valtype);
  end
end

% check that the type of the option is allowed
ok = isempty(allowedval);
for i=1:length(allowedval)
  ok = isequal(val, allowedval{i});
  if ok
    % no reason to do additional checks
    break
  end
end % for allowedtype

if ~ok
  error('the value of the option "%s" is invalid', key);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = printcell(c)
if ~isempty(c)
  s = sprintf('%s, ', c{:});
  s = sprintf('{%s}', s(1:end-2));
else
  s = '{}';
end
