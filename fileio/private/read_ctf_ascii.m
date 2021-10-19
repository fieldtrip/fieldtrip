function [file] = read_ctf_ascii(filename)

% READ_CTF_ASCII reads general data from an CTF configuration file
%
% The file should be formatted like
%    Group
%    {
%      item1 : value1a value1b value1c
%      item2 : value2a value2b value2c
%      item3 : value3a value3b value3c
%      item4 : value4a value4b value4c
%    }
%
% This fileformat structure is used in
%   params.avg
%   default.hdm
%   multiSphere.hdm
%   processing.cfg
% and maybe for other files as well.

% Copyright (C) 2003, Robert Oostenveld
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

fid = fopen_or_error(filename, 'r');

line = '';
while ischar(line)
  line = cleanline(fgetl(fid));
  if isempty(line) || (length(line)==1 && all(line==-1))
    continue
  end

  % the line is not empty, which means that we have encountered a chunck of information
  subline = cleanline(fgetl(fid));  % read the {
  subline = cleanline(fgetl(fid));  % read the first item
  while ~contains(subline, '}')
    if ~isempty(subline)
      [item, value] = strtok(subline, ':');
      value(1) = ' ';           % remove the :
      value  = strtrim(value);
      item   = strtrim(item);
      
      % the item name should be a real string, otherwise I cannot put it into the structure
      if strcmp(sprintf('%d', str2num(deblank(item))), deblank(item))
        % add something to the start of the string to distinguish it from a number
        item = ['item_' item];
      end

      % the value can be either a number or a string, and is put into the structure accordingly
      if isempty(str2num(value))
        % the value appears to be a string
        eval(sprintf('file.%s.%s = [ ''%s'' ];', line, item, value));
      else
        % the value appears to be a number or a list of numbers
        eval(sprintf('file.%s.%s = [ %s ];', line, item, value));
      end

    end
    subline = cleanline(fgetl(fid));    % read the first item
  end
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function line = cleanline(line)
  if isempty(line) || (length(line)==1 && all(line==-1))
    return
  end
  comment = strfind(line, '//');
  if ~isempty(comment)
    line(min(comment):end) = ' ';
  end
  line = strtrim(line);
