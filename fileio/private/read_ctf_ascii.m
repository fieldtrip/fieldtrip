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
% $Log: read_ctf_ascii.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.10  2008/11/14 07:49:19  roboos
% use standard matlab strtrim function instead of deblank2
%
% Revision 1.9  2006/08/29 11:08:41  roboos
% close file after reading from it
%
% Revision 1.8  2006/03/06 10:46:53  roboos
% fixed bug (again) in shortcirquited ||
%
% Revision 1.7  2006/03/06 09:40:53  roboos
% changed a | into a ||
%
% Revision 1.6  2006/02/24 15:43:43  roboos
% fixed bug in shortcirquited || which was not boolean on both sides
%
% Revision 1.5  2006/02/09 08:36:33  roboos
% changed single | into || since it is a boolean evaluation (thanks to Tom)
%
% Revision 1.4  2004/08/02 13:06:09  roboos
% fixed bug that occurred for hdm files based on markers: structs cannot have fields that have a name which consists of a single number, so in that case use the fieldname "item_1" etc.
%
% Revision 1.3  2004/06/28 07:32:28  roberto
% added warning off/on around text reading and conversion
%
% Revision 1.2  2003/04/17 12:38:08  roberto
% *** empty log message ***
%
% Revision 1.1  2003/03/24 12:30:42  roberto
% new implementation
%

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

