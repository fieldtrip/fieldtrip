function [val] = filetype_check_extension(filename, ext)

% FILETYPE_CHECK_EXTENSION helper function to determine the file type
% by performing a case insensitive string comparison of the extension.

% Copyright (C) 2003-2006 Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

