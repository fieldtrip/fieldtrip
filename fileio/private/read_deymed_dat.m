function dat = read_deymed_dat(filename, hdr, begsample, endsample)

% READ_DEYMED_DAT reads EEG data from the Deymed Truescan file format
%
% Use as
%   dat = read_deymed_dat(filename, hdr, begsample, endsample)
%
% See also READ_DEYMED_INI

% Copyright (C) 2013, Robert Oostenveld
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

[p, f, x] = fileparts(filename);
headerfile = fullfile(p, [f '.ini']);
if ~exist(headerfile, 'file')
  headerfile = fullfile(p, [f '.Ini']);
end
datafile = fullfile(p, [f '.154576dat']);
if ~exist(datafile, 'file')
  datafile = fullfile(p, [f '.Dat']);
end

if nargin<2 || isempty(hdr)
  hdr = read_deymed_ini(headerfile);
end
if nargin<3 || isempty(begsample)
  begsample = 1;
end
if nargin<4 || isempty(endsample)
  endsample = hdr.nSamples;
end

fid = fopen_or_error(datafile, 'rb', 'ieee-le');
offset = 512;                                 % for the general header
offset = offset + hdr.nChans*1024;            % for the channel headers
offset = offset + (begsample-1)*2*hdr.nChans; % for the begin sample
fseek(fid, offset, 'cof');
dat = fread(fid, [hdr.nChans (endsample-begsample+1)], 'int16');
fclose(fid);

% calibrate to convert integers to uV
dat = dat*hdr.calib;



