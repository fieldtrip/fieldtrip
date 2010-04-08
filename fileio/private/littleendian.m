function val = littleendian;

% LITTLEENDIAN returns 1 (true) on a little endian machine, e.g. with an
% Intel or AMD, or 0 (false) otherwise
%
% Example
%   if (littleendian)
%     % do something, e.g. swap some bytes
%    end
%
% See also BIGENDIAN, SWAPBYTES, TYPECAST

% Copyrigth (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

val = (typecast(uint8([0 1]), 'uint16')==256);
