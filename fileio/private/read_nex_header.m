function [hdr] = read_nex_header(filename)

% READ_NEX_HEADER for Plexon *.nex file
%
% Use as
%   [hdr] = read_nex_header(filename)
%
% See also RAD_NEX_DATA, READ_NEX_EVENT

% Copyright (C) 2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

fid = fopen(filename, 'r', 'ieee-le');

% reading the file header
filheader.magicnumber     = fread(fid,4,'uint8=>char')';
filheader.version         = fread(fid,1,'int32');
filheader.comment         = fread(fid,256,'uint8=>char')';
filheader.frequency       = fread(fid,1,'double');
filheader.begvar          = fread(fid,1,'int32');
filheader.endvar          = fread(fid,1,'int32');
filheader.numvar          = fread(fid,1,'int32');
filheader.nextfileheader  = fread(fid,1,'int32');
filheader.padding         = fread(fid,256,'uint8=>char')';

% reading the variable headers
for varlop=1:filheader.numvar
  varheader(varlop).typ        = fread(fid,1,'int32');
  varheader(varlop).version    = fread(fid,1,'int32');
  varheader(varlop).nam        = fread(fid,64,'uint8=>char')';
  varheader(varlop).offset     = fread(fid,1,'int32');
  varheader(varlop).cnt        = fread(fid,1,'int32');
  varheader(varlop).wirenumber = fread(fid,1,'int32');
  varheader(varlop).unitnumber = fread(fid,1,'int32');
  varheader(varlop).gain       = fread(fid,1,'int32');
  varheader(varlop).filter     = fread(fid,1,'int32');
  varheader(varlop).xpos       = fread(fid,1,'double');
  varheader(varlop).ypos       = fread(fid,1,'double');
  varheader(varlop).wfrequency = fread(fid,1,'double');
  varheader(varlop).adtomv     = fread(fid,1,'double');
  varheader(varlop).numsmp     = fread(fid,1,'int32');
  varheader(varlop).nummrk     = fread(fid,1,'int32');
  varheader(varlop).mrklen     = fread(fid,1,'int32');
  padding                      = fread(fid,68,'uint8=>char')';
end

status = fclose(fid);

% put them together into one struct
hdr.fil       = filename;
hdr.filheader = filheader;
hdr.varheader = varheader;
