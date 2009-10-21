function [val] = filetype_check_extension(filename, ext)

% FILETYPE_CHECK_EXTENSION helper function to determine the file type
% by performing a case insensitive string comparison of the extension.

% Copyright (C) 2003-2006 Robert Oostenveld
%
% $Log: filetype_check_extension.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.3  2008/07/01 13:35:47  roboos
% look at number of elements, not length
%
% Revision 1.2  2007/03/21 17:22:08  roboos
% use recursion instead of copying the same code twice
%
% Revision 1.1  2006/08/28 10:10:11  roboos
% moved subfunctions for filetype detection into seperate files
%

if iscell(filename)
  % compare the extension of multiple files
  val = zeros(size(filename));
  for i=1:length(filename)
    val(i) = filetype_check_extension(filename{i}, ext);
  end
else
  % compare the extension of a single file
  if numel(filename)<numel(ext)
    val = 0;
  else
    val = strcmpi(filename((end-length(ext)+1):end), ext);
  end
end
return

