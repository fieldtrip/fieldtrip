function [val] = read_asa(filename, elem, format, number, token)

% READ_ASA reads a specified element from an ASA file
%
% val = read_asa(filename, element, type, number)
%
% where the element is a string such as
%	NumberSlices
%	NumberPositions
%	Rows
%	Columns
%	etc.
%
% and format specifies the datatype according to
%	%d	(integer value)
%	%f	(floating point value)
%	%s	(string)
%
% number is optional to specify how many lines of data should be read
% The default is 1 for strings and Inf for numbers.
%
% token is optional to specifiy a character that separates the values from
% anything not wanted.

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_asa.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.7  2008/11/17 17:38:43  roboos
% replaced line==-1 with isequal, to work around && bug
%
% Revision 1.6  2008/11/14 07:49:19  roboos
% use standard matlab strtrim function instead of deblank2
%
% Revision 1.5  2008/11/14 07:21:45  roboos
% newer ASA versions write the labels in front of the positions, like "FPz: 10.4 1.3 -2"
% added support for this, thanks to Thomas Hartmann
% use strcmpi instead of strcmp(lower())
%
% Revision 1.4  2005/11/30 11:38:33  roboos
% fixed bug: close file before returning to calling function (thanks to Gijs)
%
% Revision 1.3  2004/03/29 15:15:12  roberto
% unknown change, seems related to handling of empty lines ?
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

fid = fopen(filename, 'rt');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

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
