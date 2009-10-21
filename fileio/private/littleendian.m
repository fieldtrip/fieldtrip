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
% $Log: littleendian.m,v $
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.1  2007/01/04 12:10:42  roboos
% new implementation
%

val = (typecast(uint8([0 1]), 'uint16')==256);
