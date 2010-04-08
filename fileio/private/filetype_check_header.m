function [val] = filetype_check_header(filename, head, offset)

% FILETYPE_CHECK_HEADER helper function to determine the file type
% by reading the first number of bytes of a file and comparing them
% to a known string (c.f. magic number).

% Copyright (C) 2003-2006 Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if nargin<3
  offset = 0;
end

if iscell(filename)
  % compare the header of multiple files
  val = zeros(size(filename));
  for i=1:length(filename)
    val(i) = filetype_check_header(filename{i}, head, offset);
  end
elseif isdir(filename)
  % a directory cannot have a header
  val = 0;
else
  % read the first few bytes from the file and compare them to the desired header
  fid = fopen(filename, 'rb');
  if fid<0
    warning(sprintf('could not open %s', filename));
    val = 0;
  else
    fseek(fid, offset, 'cof');
    if iscell(head)
      for i=1:length(head)
        len(i) = length(head{i});
      end
      [str, siz] = fread(fid, max(len), 'uint8=>char');
      fclose(fid);
      for i=1:length(head)
        val = strncmp(str, head{i}, len(i));
        if val
          break
        end
      end
    else
      [str, siz] = fread(fid, length(head), 'uint8=>char');
      fclose(fid);
      if siz~=length(head)
        warning(sprintf('could not read the header from %s', filename));
        val = 0;
      else
        val = all(str(:)==head(:));
      end
    end
  end
end
return

