function [file] = read_bti_ascii(filename)

% READ_BTI_ASCII reads general data from a BTI configuration file
%
% The file should be formatted like
%    Group: 
%      item1 : value1a value1b value1c 
%      item2 : value2a value2b value2c 
%      item3 : value3a value3b value3c 
%      item4 : value4a value4b value4c 
%

% Copyright (C) 2004, Robert Oostenveld
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

fid = fopen(filename, 'r');
if fid==-1
  ft_error(sprintf('could not open file %s', filename));
end

line = '';
while ischar(line)
  line = cleanline(fgetl(fid));
  
  if isempty(line) || line==-1 || isempty(findstr(line, ':'))
    continue
  end
  
  % the line is not empty, which means that we have encountered a chunck of information
  if findstr(line, ':')~=length(line)
    [item, value] = strtok(line, ':');
    value(1) = ' ';         % remove the :
    value  = strtrim(value);
    item   = strtrim(item);
    item(findstr(item, '.')) = '_';
    item(findstr(item, ' ')) = '_';
    if ischar(item)
      eval(sprintf('file.%s = ''%s'';', item, value));
    else
      eval(sprintf('file.%s = %s;', item, value));
    end  
  else
    subline = cleanline(fgetl(fid));
    error, the rest has not been implemented (yet)
    
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function line = cleanline(line)
if isempty(line) || line==-1
  return
end
comment = findstr(line, '//');
if ~isempty(comment)
  line(min(comment):end) = ' ';
end
line = strtrim(line);

