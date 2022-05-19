function ft_write_json(filename, json, varargin)

% FT_WRITE_JSON writes a MATLAB structure to a JSON file. Compared to the builtin
% MATLAB function, this implementation deals a bit different with missing values,
% booleans, and NaNs, and results in a more human-readable file.
%
% Use as
%   ft_write_json(filename, struct)
%
% See also FT_READ_JSON, JSONDECODE, JSONENCODE

% Copyright (C) 2018-2021, Robert Oostenveld
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

% get the options
sort = ft_getopt(varargin, 'sort', true); % can be true/false, or 'yes'/'no'

ft_info('writing ''%s''\n', filename);
json = remove_empty(json);
if istrue(sort)
 json = sort_fields(json);
end
json = ft_struct2char(json); % convert strings into char-arrays
ft_hastoolbox('jsonlab', 1);
% write nan as 'n/a'
% write boolean as True/False
str = savejson('', json, 'NaN', '"n/a"', 'ParseLogical', true);
fid = fopen_or_error(filename, 'w');
fwrite(fid, str);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this is shared with DATA2BIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = remove_empty(s)
if isempty(s)
  return
elseif isstruct(s)
  fn = fieldnames(s);
  fn = fn(structfun(@isempty, s));
  s = removefields(s, fn);
elseif istable(s)
  remove = false(1,size(s,2));
  for i=1:size(s,2)
    % find columns that are non-numeric and where all elements are []
    remove(i) = ~isnumeric(s{:,i}) && all(cellfun(@isempty, s{:,i}));
  end
  s = s(:,~remove);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this is shared with DATA2BIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = sort_fields(x)
fn = fieldnames(x);
fn = sort(fn);
y = struct();
for i=1:numel(fn)
  y.(fn{i}) = x.(fn{i});
end
