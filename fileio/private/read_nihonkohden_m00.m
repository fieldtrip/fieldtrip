function [dat] = read_nihonkohden_m00(filename, begsample, endsample)

% READ_NIHONKOHDEN_M00 reads data Nihon Kohden EEG system and converts the
% floating point ASCII values to 16-bit signed integers from *.m00 file
% format. This function is an adaptation of convert_nkascii2mat.m written
% by Timothy Ellmore
%
% Use as
%   dat = read_nihonkohden_m00(filename)
%
% This returns a header structure with
%   hdr.nSamples      number of data points per channel
%   hdr.nchannels     number of channels sampled during recording
%   hdr.bsweep        begin sweep (ms)
%   hdr.sampintms     sampling interval in ms
%   hdr.binsuV        number of bins per microvolts
%   hdr.start_time    starting time of recording
%   hdr.label         char array of channel names
%
% See also READ_NIHONKOHDEN_HDR

% Copyright (C) 2016, Diego Lozano-Soldevilla
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

hdr = read_nihonkohden_hdr(filename);
mult = 10;% a multiplier; divide stored value by this to floats
nsmp = (endsample-begsample+1);

% the name of the ascii file
fid=fopen(filename,'r');

% skip header information
hline1 = fgets(fid);
hline2 = fgets(fid);

fprintf('\nWill read from %d to %d sample points: %d seconds \n',begsample,endsample,nsmp/hdr.Fs);

dc = textscan(fid,'%f',hdr.nChans*nsmp,'delimiter','\n','headerlines',begsample-1);
df = dc{1};clear dc;
df = df.*mult; % multiply by 10, don't forget to divide when reading in data!

fprintf('converting ASCII floats to 16-bit signed ints');

dat = reshape(df,hdr.nChans,nsmp);
clear df;
dat = int16(dat)./mult;
dat = double(dat);
% should we double it? No clue why the int16 and the multiplier.
% Later convert data has to be multiplied to transform it into microvolts
% nkdata.eeg=nkdata.eeg./nkdata.multiplier; (Diego)
fclose(fid);
