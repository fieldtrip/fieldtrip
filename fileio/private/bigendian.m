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
% $Log: bigendian.m,v $
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.1  2007/01/04 12:10:42  roboos
% new implementation
%

val = (typecast(uint8([0 1]), 'uint16')==1);
