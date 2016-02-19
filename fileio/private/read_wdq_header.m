function [hdr] = read_wdq_header(filename)

% READ_WDQ_HEADER reads header information from wdq files
%
% Use as
%  [hdr] = read_wdq_header(filename)
%

% Copyright (C) 2011, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% information about how to interpret the file are taken from the document
% 'CODAS data storage format'

fid = fopen(filename, 'r');
%tmp = fread(fid, 110); % the first 110 bytes are always interpretable equally
hdr = [];
hdr.nchan         = fread(fid, 1, 'uint16=>uint16');
% hdr.nchan         = hdr.nchan - bitand(hdr.nchan, 32);
hdr.nchan         = double(bitand(hdr.nchan, 255)); % only the lowest bits seem to be coding the number of channels
i2                = fread(fid, 1, 'uint16');
hdr.offset        = fread(fid, 1, 'uint8');
hdr.nbyteschanhdr = fread(fid, 1, 'uint8');
hdr.nbyteshdr     = fread(fid, 1, 'int16');
hdr.nbytesdat     = fread(fid, 1, 'uint32'); %FIXME only for unpacked files
hdr.nbytestrailer1 = fread(fid, 1, 'uint32');
hdr.nbytestrailer2 = fread(fid, 1, 'uint16');
i9  = fread(fid, 1, 'uint16'); % not used
i10 = fread(fid, 1, 'uint16'); % not used
i11 = fread(fid, 1, 'int16');
i12 = fread(fid, 4, 'int8');
i13 = fread(fid, 1, 'double');

% FIXME there seems to be a discrepancy here, according to the
% documentation it should be the next (commented) line. according to the
% file I have at hand it should be the uncommented line.
%hdr.fsample        = double(hdr.nchan)./i13;
hdr.fsample = 1./i13;
i14 = fread(fid, 1, 'uint32');
i15 = fread(fid, 1, 'uint32');
i16 = fread(fid, 1, 'uint32');
i17 = fread(fid, 1, 'uint32');
i18 = fread(fid, 1, 'uint32');
i19 = fread(fid, 2, 'uint16');
i20 = fread(fid, 1, 'int16');
i21 = fread(fid, 1, 'int16');
i22 = fread(fid, 1, 'int8');
i23 = fread(fid, 1, 'int8');
i24 = fread(fid, 1, 'int8');
i25 = fread(fid, 1, 'int8');
i26 = fread(fid, 32, 'int8');
i27 = fread(fid, 1, 'uint16');
i28 = fread(fid, 1, 'uint8');
empty = fread(fid, 1, 'uint8');
i29 = fread(fid, 1, 'int8');
i30 = fread(fid, 1, 'int8');
i31 = fread(fid, 1, 'int8');
empty = fread(fid, 1, 'uint8');
i32 = fread(fid, 1, 'int8');
i33 = fread(fid, 1, 'int8');

% interpret relevant stuff from i27
% other stuff may or may not be relevant
hdr.ispacked  = bitand(i27, 2^15)>0;
hdr.lowchan   = bitand(i27, 2^9) >0;

% get channel info
for k = 1:(hdr.nbyteshdr-hdr.offset-2)/hdr.nbyteschanhdr
  fseek(fid, hdr.offset + (k-1)*hdr.nbyteschanhdr, 'bof');
  chanhdr(k).scalingslope     = fread(fid, 1, 'float');
  chanhdr(k).scalingintercept = fread(fid, 1, 'float');
  chanhdr(k).scale            = fread(fid, 1, 'double');
  chanhdr(k).intercept        = fread(fid, 1, 'double');
  chanhdr(k).unit             = fread(fid, 1, 'uint8=>char');
  fread(fid, 1, 'uint8');
  fread(fid, 1, 'uint8'); % may be relevant for packed files
  chanhdr(k).channr = fread(fid, 1, 'uint8');
  chanhdr(k).label  = num2str(k);
end
hdr.chanhdr = chanhdr;

fseek(fid, hdr.nbyteshdr + hdr.nbytesdat, 'bof');
x = fread(fid, hdr.nbytestrailer1, 'uint8');
hdr.x = x; % what shall we do with this?

fseek(fid, hdr.nbyteshdr + hdr.nbytesdat + hdr.nbytestrailer1, 'bof');
annot = fread(fid, hdr.nbytestrailer2, 'uint8');
sel   = find(annot==0);
indx  = find(annot>0, 1, 'first');
sel   = [0;sel];
for k = 1:numel(sel)-1
  hdr.chanhdr(k).annot = char(annot(sel(k)+1:sel(k+1)-1))';  
end

fclose(fid);
