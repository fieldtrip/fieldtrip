function [dat] = read_plexon_ddt(filename, begsample, endsample)

% READ_PLEXON_DDT reads header or data from a Plexon *.ddt file,
% which is a Plexon continuous data file optimized for continuous
% (streaming) recording where every channel is continuously recorded
% without gaps and the recording includes any dead time between spikes.
% 
% Use as
%   [hdr] = read_plexon_ddt(filename)
%   [dat] = read_plexon_ddt(filename, begsample, endsample)
%
% samples start counting at 1
% returned values are in mV

% Copyright (C) 2005-2007, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 100: Samples are assumed to be 12 bits.  All channels have the same NIDAQ gain, and preamp 
% gain is not saved. 
%     int     Version;         // =100 
%     int     DataOffset;      // Offset into the file where the data starts 
%     double  Freq;            // Digitization frequency 
%     int     NChannels;       // Number of channels 
%     int     Year;            // Time/date when the data was acquired 
%     int     Month; 
%     int     Day; 
%     int     Hour; 
%     int     Minute; 
%     int     Second; 
%     int     Gain;            // NIDAQ gain (the same gain for all channels 
%     char    Comment[128];    // User-supplied comment  
%     unsigned char Padding[256]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat = [];
fid = fopen(filename, 'rb', 'ieee-le');

Version    = fread(fid, 1, 'int'    ); 
dat.Version    = Version;

% this is common for all subformats
dat.DataOffset = fread(fid, 1, 'int'    );
dat.Freq       = fread(fid, 1, 'double' );
dat.NChannels  = fread(fid, 1, 'int'    );
dat.Year       = fread(fid, 1, 'int'    );
dat.Month      = fread(fid, 1, 'int'    );
dat.Day        = fread(fid, 1, 'int'    );
dat.Hour       = fread(fid, 1, 'int'    );
dat.Minute     = fread(fid, 1, 'int'    );
dat.Second     = fread(fid, 1, 'int'    );
dat.Gain       = fread(fid, 1, 'int'    );
dat.Comment        = char(fread(fid, 128, 'char')');

if Version==100
  % the remainder of the header contains bytes for padding
  dat.Padding    = fread(fid, 256, 'char' );
elseif Version==101
  % Version 101: A field was added to indicate bits-per-sample (12 or 16).  All channels have the same NIDAQ
  % gain, and preamp gain is not saved.
  dat.BitsPerSample = fread(fid,   1, 'char');
  dat.Padding       = fread(fid, 255, 'char');
elseif Version==102
  % Version 102: A byte array of per-channel NIDAQ gains was added.  The Gain field now indicates the
  % preamp gain,   a value of 1 usually indicating that no preamp gain was specified in the application that
  % recorded the DDT.
  dat.BitsPerSample = fread(fid,   1, 'char');
  dat.ChannelGain   = fread(fid,  64, 'char');
  dat.Padding       = fread(fid, 191, 'char');
elseif Version==103
  % Version 103: A field was added to specify ADC maximal input voltage. The value is integer number in
  % millivolts.
  dat.BitsPerSample  = fread(fid,   1, 'char');
  dat.ChannelGain    = fread(fid,  64, 'char');  
  dat.MaxMagnitudeMV = fread(fid,   1, 'short');
  dat.Padding        = fread(fid, 189, 'char');
else
  error('unsupported version of ddt file');
end

% determine the number of samples by looking at the length of the datafile
offset = ftell(fid);
fseek(fid, 0, 'eof');
dat.NSamples = (ftell(fid) - offset) / (2 * dat.NChannels);

if nargin==1
  % do not read any data, only return the header
  fclose(fid);
  return
end

% A DDT file contains a file header followed by an array of A/D samples stored as 16-bit integers regardless 
% of whether values are 12 or 16-bits. The samples are multiplexed.

% read the desired samples of the data
offset = dat.DataOffset + (begsample-1)*dat.NChannels;
nsampl = endsample - begsample + 1;
fseek(fid, offset, 'bof');
dat.data = fread(fid, [dat.NChannels nsampl], 'int16');

fclose(fid);

% apply the calibration to obtain the signal in physical units
if Version==100
  PreAmpGain = 1000;
  dat.data = dat.data * 5000 / (2048 * dat.Gain * PreAmpGain);
elseif Version==101
  PreAmpGain = 1000;
  dat.data = dat.data * 5000 / (0.5 * (2^dat.BitsPerSample) * dat.Gain * PreAmpGain);
elseif Version==102
  for i=1:dat.NChannels
    dat.data(i,:) = dat.data(i,:) * 5000 ./ (0.5 * (2^dat.BitsPerSample) * dat.ChannelGain(i));
  end
elseif Version==103
  % I am not sure whether the calibration like this is correct, since the
  % Plexon documentation does not explicitely specify how to do it for the
  % 104 file format. I presume that it is identical to the 103 format.
  for i=1:dat.NChannels
    dat.data(i,:) = dat.data(i,:) * 5000 ./ (0.5 * (2^dat.BitsPerSample) * dat.ChannelGain(i));
  end
end

