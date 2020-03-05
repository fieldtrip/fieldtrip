function fid = fopen_or_error(filename, mode, varargin)

% FOPEN_OR_ERROR Opens a file, like fopen, but throws an exception if the open failed.
%
% This keeps you from having to write "if fid < 0; error(...)" everywhere
% you do an fopen.
%
% See also ISDIR_OR_MKDIR

if nargin < 2 || isempty(mode)
  mode = 'r';
end

[fid, msg] = fopen(filename, mode, varargin{:});

if fid < 0
  if any(mode == 'w')
    mode_descr = 'writing';
  elseif any(mode == 'a')
    mode_descr = 'appending';
  else
    mode_descr = 'reading';
  end  
  ft_error(sprintf('could not open file %s for %s: %s', ...
    filename, mode_descr, msg));
end

