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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
if isempty(line) | line==-1
  return
end
comment = findstr(line, '//');
if ~isempty(comment)
  line(min(comment):end) = ' ';
end
line = strtrim(line);

