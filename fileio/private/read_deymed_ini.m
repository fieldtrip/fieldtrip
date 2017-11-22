function hdr = read_deymed_ini(filename)

% READ_DEYMED_INI reads EEG data from the Deymed Truescan file format
%
% Use as
%   hdr = read_deymed_ini(filename)
%
% See also READ_DEYMED_DAT

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

fid = fopen(datafile, 'rb');
hdr.id        = fread(fid, [1, 12], 'uint8=>char');
hdr.lastname  = fread(fid, [1, 12], 'uint8=>char');
hdr.firstname = fread(fid, [1, 12], 'uint8=>char');
hdr.date      = fread(fid, [1, 12], 'uint8=>char');
hdr.time      = fread(fid, [1, 12], 'uint8=>char');
hdr.Fs        = fread(fid, [1, 1],  'int16');
hdr.nChans    = fread(fid, [1, 1],  'int8');

dummy         = fread(fid, [1, 1],  'int8'); % this is 4 in the example file
if dummy~=4
  ft_warning('deviant header format detected');
end

% A numerical comparison between a deymed dat file and a BESA-exported edf version of the
% same file showed a calibration difference of 1/4. So this might be a calibration value.
hdr.calib = 1/dummy;

hdr.label = {};
for i=1:hdr.nChans
  str = fread(fid, [1, 6], 'uint8=>char');
  str(str==0)  = ' ';
  hdr.label{i} = strtrim(str);
end
fclose(fid);

% determine the number of samples
d = dir(datafile);
hdr.nSamples = (d.bytes-512-1024*hdr.nChans)/(2*hdr.nChans);

% also store some details from the ini file
try
  hdr.ini.patient.id        = inifile(headerfile, 'read', {'Pacient', '', 'ID'});
  hdr.ini.patient.lastname  = inifile(headerfile, 'read', {'Pacient', '', 'Prijmeni'});
  hdr.ini.patient.firstname = inifile(headerfile, 'read', {'Pacient', '', 'Jmeno'});
  hdr.ini.record.user       = inifile(headerfile, 'read', {'Record',  '', 'AcqUser'});
  hdr.ini.record.type       = inifile(headerfile, 'read', {'Record',  '', 'Type'});
end
