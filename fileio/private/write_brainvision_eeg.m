function write_brainvision_eeg(filename, hdr, dat, event)

% WRITE_BRAINVISION_EEG exports continuous EEG data to a BrainVision *.eeg
% and corresponding *.vhdr file. The samples in the exported file are
% multiplexed and stored in ieee-le float32 format.
%
% Use as
%   write_brainvision_eeg(filename, hdr, dat, evt)
%
% See also READ_BRAINVISION_EEG, READ_BRAINVISION_VHDR, READ_BRAINVISION_VMRK

% Copyright (C) 2007-2014, Robert Oostenveld
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

if nargin<4
  event = [];
end

if length(size(dat))>2
  ntrl  = size(dat,1);
  nchan = size(dat,2);
  nsmp  = size(dat,3);
else
  nchan = size(dat,1);
  nsmp  = size(dat,2);
end

if hdr.nChans~=nchan
  ft_error('number of channels in in header does not match with the data');
end

% this is the only supported data format
hdr.DataFormat      = 'BINARY';
hdr.DataOrientation = 'MULTIPLEXED';
hdr.BinaryFormat    = 'IEEE_FLOAT_32';
hdr.resolution      = ones(size(hdr.label));  % no additional calibration needed, since float32

% determine the filenames
[p, f, x] = fileparts(filename);
headerfile = fullfile(p, [f '.vhdr']);
markerfile = fullfile(p, [f '.vmrk']);
datafile   = fullfile(p, [f '.eeg']);

% the files internally refer to each other, this should be without the leading directory
headerfile_without_path = [f '.vhdr'];
markerfile_without_path = [f '.vmrk'];
datafile_without_path   = [f '.eeg'];


% open the data file and write the binary data
fid = fopen_or_error(datafile, 'wb', 'ieee-le');
if length(size(dat))>2
  ft_warning('writing segmented data as if it were continuous');
  for i=1:ntrl
    fwrite(fid, squeeze(dat(i,:,:)), 'float32');
  end
else
  fwrite(fid, dat, 'float32');
end

fclose(fid);

% open the header file and write the ascii header information
fid = fopen_or_error(headerfile, 'wb');
fprintf(fid, 'Brain Vision Data Exchange Header File Version 1.0\r\n');
fprintf(fid, '; Data created by FieldTrip\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[Common Infos]\r\n');
fprintf(fid, 'DataFile=%s\r\n',          datafile_without_path);
if ~isempty(markerfile)
  fprintf(fid, 'MarkerFile=%s\r\n',      markerfile_without_path);
end
fprintf(fid, 'DataFormat=%s\r\n',        hdr.DataFormat);
fprintf(fid, 'DataOrientation=%s\r\n',   hdr.DataOrientation);
fprintf(fid, 'NumberOfChannels=%d\r\n',  hdr.nChans);
% Sampling interval in microseconds
fprintf(fid, 'SamplingInterval=%d\r\n',  1e6/hdr.Fs);
fprintf(fid, '\r\n');
fprintf(fid, '[Binary Infos]\r\n');
fprintf(fid, 'BinaryFormat=%s\r\n',      hdr.BinaryFormat);
fprintf(fid, '\r\n');
fprintf(fid, '[Channel Infos]\r\n');
% Each entry: Ch<Channel number>=<Name>,<Reference channel name>,<Resolution in microvolts>,<Future extensions>...
% Fields are delimited by commas, some fields might be omitted (empty).
% Commas in channel names should be coded as "\1", but are not supported here
for i=1:hdr.nChans
  fprintf(fid, 'Ch%d=%s,,%g\r\n', i, hdr.label{i}, hdr.resolution(i));
end
fclose(fid);

% open the marker file and write the ascii marker information
fid = fopen_or_error(markerfile, 'wb');
fprintf(fid, 'Brain Vision Data Exchange Marker File, Version 1.0\n');
fprintf(fid, '\n');
fprintf(fid, '[Common Infos]\n');
fprintf(fid, 'Codepage=UTF-8\n');
fprintf(fid, 'DataFile=%s\n', datafile_without_path);
fprintf(fid, '\n');
fprintf(fid, '[Marker Infos]\n');
fprintf(fid, '; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,\n');
fprintf(fid, '; <Size in data points>, <Channel number (0 = marker is related to all channels)>\n');
fprintf(fid, '; Fields are delimited by commas, some fields might be omitted (empty).\n');
fprintf(fid, '; Commas in type or description text are coded as "\1".\n');
for i=1:length(event)
  type  = event(i).type;          % type is always a string
  descr = event(i).value;         % value can be empty, string or numeric
  if isempty(descr)
    descr = '';
  elseif isnumeric(descr)
    descr = num2str(descr);
  end
  pos = num2str(event(i).sample); % sample is always numeric, hence convert to string
  siz = event(i).duration;        % duration can be empty or numeric
  if isempty(siz)
    siz = '1';
  else
    siz = num2str(siz);
  end
  chan = '0';
  fprintf(fid, 'Mk%d=%s,%s,%s,%s,%s\n', i, type, descr, pos, siz, chan);
end
fclose(fid);
