function fid = fopen_or_error(filename, permission, varargin)

% FOPEN_OR_ERROR Opens a file, like fopen, but throws an exception if the open failed.
%
% This keeps you from having to write "if fid < 0; error(...)" everywhere
% you do an fopen.
%
% See also FOPEN, ISDIR_OR_MKDIR

if nargin < 2 || isempty(permission)
  permission = 'r';
end

[fid, msg] = fopen(filename, permission, varargin{:});

if fid < 0
  if any(permission == 'w')
    mode_descr = 'writing';
  elseif any(permission == 'a')
    mode_descr = 'appending';
  else
    mode_descr = 'reading';
  end  
  ft_error(sprintf('could not open file %s for %s: %s', ...
    filename, mode_descr, msg));
end

