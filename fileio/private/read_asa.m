function [val] = read_asa(filename, elem, format, number, token)

% READ_ASA reads a specified element from an ASA file
%
% val = read_asa(filename, element, type, number)
%
% where the element is a string such as
%   NumberSlices
%   NumberPositions
%   Rows
%   Columns
%   etc.
%
% and format specifies the datatype according to
%   %d  (integer value)
%   %f  (floating point value)
%   %s  (string)
%
% number is optional to specify how many lines of data should be read
% The default is 1 for strings and Inf for numbers.
%
% token is optional to specifiy a character that separates the values from
% anything not wanted.

% Copyright (C) 2002-2012, Robert Oostenveld
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

fid = fopen_or_error(filename, 'rt');

if nargin<4
  if strcmp(format, '%s')
    number = 1;
  else
    number = Inf;
  end
end

if nargin<5
  token = '';
end


val = [];
elem = strtrim(lower(elem));

while (1)
  line = fgetl(fid);
  if ~isempty(line) && isequal(line, -1)
    % prematurely reached end of file
    fclose(fid);
    return
  end
  line = strtrim(line);
  lower_line = lower(line);
  if strmatch(elem, lower_line)
    data = line((length(elem)+1):end);
    break
  end
end

while isempty(data)
  line = fgetl(fid);
  if isequal(line, -1)
    % prematurely reached end of file
    fclose(fid);
    return
  end
  data = strtrim(line);
end

if strcmp(format, '%s')
  if number==1
    % interpret the data as a single string, create char-array
    val = detoken(strtrim(data), token);
    if val(1)=='='
      val = val(2:end); % remove the trailing =
    end
    fclose(fid);
    return
  end
  % interpret the data as a single string, create cell-array
  val{1} = detoken(strtrim(data), token);
  count = 1;
  % read the remaining strings
  while count<number
    line = fgetl(fid);
    if ~isempty(line) && isequal(line, -1)
      fclose(fid);
      return
    end
    tmp = sscanf(line, format);
    if isempty(tmp)
      fclose(fid);
      return
    else
      count = count + 1;
      val{count} = detoken(strtrim(line), token);
    end
  end

else
  % interpret the data as numeric, create numeric array
  count = 1;
  data = sscanf(detoken(data, token), format)';
  if isempty(data),
    fclose(fid);
    return
  else
    val(count,:) = data;
  end
  % read remaining numeric data
  while count<number
    line = fgetl(fid);
    if ~isempty(line) && isequal(line, -1)
      fclose(fid);
      return
    end
    data = sscanf(detoken(line, token), format)';
    if isempty(data)
      fclose(fid);
      return
    else
      count = count+1;
      val(count,:) = data;
    end
  end
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = detoken(in, token)
if isempty(token)
  out = in;
  return;
end

[tok rem] = strtok(in, token);
if isempty(rem)
  out = in;
  return;
else
  out = strtok(rem, token);
  return
end
