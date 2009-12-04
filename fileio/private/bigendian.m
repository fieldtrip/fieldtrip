function val = bigendian;

% BIGENDIAN returns 1 (true) on a big endian machine, e.g. with a SUN Sparc
% or Apple G4 processor, or 0 (false) otherwise
%
% Example
%   if (bigendian)
%     % do something, e.g. swap some bytes
%    end
%
% See also LITTLEENDIAN, SWAPBYTES, TYPECAST

% Copyrigth (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

val = (typecast(uint8([0 1]), 'uint16')==1);
