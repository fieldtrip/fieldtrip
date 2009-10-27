function [file] = read_ctf_ascii(filename);

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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fid = fopen(filename, 'r');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

line = '';
while ischar(line)
  line = cleanline(fgetl(fid));
  if isempty(line) || (length(line)==1 && all(line==-1))
    continue
  end

  % the line is not empty, which means that we have encountered a chunck of information
  subline = cleanline(fgetl(fid));	% read the {
  subline = cleanline(fgetl(fid));	% read the first item
  while isempty(findstr(subline, '}'))
    if ~isempty(subline)
      [item, value] = strtok(subline, ':');
      value(1) = ' ';			% remove the :
      value  = strtrim(value);
      item   = strtrim(item);
      warning off
      
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
      warning on
    end
    subline = cleanline(fgetl(fid));	% read the first item
  end
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function line = cleanline(line)
  if isempty(line) || (length(line)==1 && all(line==-1))
    return
  end
  comment = findstr(line, '//');
  if ~isempty(comment)
    line(min(comment):end) = ' ';
  end
  line = strtrim(line);

