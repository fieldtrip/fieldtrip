function [file] = read_bti_ascii(filename);

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
% $Log: read_bti_ascii.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.2  2008/11/14 07:49:19  roboos
% use standard matlab strtrim function instead of deblank2
%
% Revision 1.1  2005/04/18 13:47:05  roboos
% added some old and infrequently used functions to the cvs repository
%
% Revision 1.1  2005/04/18 13:43:34  roboos
% included some old functions in the cvs repository, this ensures consistency of the functions between the different network locations
%

fid = fopen(filename, 'r');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

line = '';
while ischar(line)
  line = cleanline(fgetl(fid))
  
  if isempty(line) | line==-1 | isempty(findstr(line, ':'))
    continue
  end
  
  % the line is not empty, which means that we have encountered a chunck of information
  if findstr(line, ':')~=length(line)
    [item, value] = strtok(line, ':');
    value(1) = ' ';			% remove the :
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
if isempty(line) | line==-1
  return
end
comment = findstr(line, '//');
if ~isempty(comment)
  line(min(comment):end) = ' ';
end
line = strtrim(line);

