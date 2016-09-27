function [hdr] = read_nihonkohden_hdr(filename)

% READ_NIHONKOHDEN_HDR reads the known items from the Nihon Kohden EEG
% system header file (*.m00) and returns them in a structure. This reading
% function is an adaptation of get_nkheader.m written by Timothy Ellmore
%
% Use as
%   hdr = read_nihonkoden_hdr(filename)
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
% See also READ_NIHONKOHDEN_M00

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

% the name of the ascii file
fid=fopen(filename,'r');

%header line 1: acquisition information
hline1 = fgets(fid);

% get the six different fields in the header, store in cell array
%hd=textscan(hline1,'%s %s %s %s %s %s');
[hd]=strread(hline1,'%s');

fprintf('\nNihon Kohden ASCII EEG header fields:');
fprintf('\n-------------------------------------');

% number of timepoints
[txt,hdr.nSamples]=strread(char(hd{1}),'%s%d','delimiter','=');
fprintf('\n%s is %d',char(txt),hdr.nSamples);

% number of channels sampled
[txt,hdr.nChans]=strread(char(hd{2}),'%s%d','delimiter','=');
fprintf('\nNumber of %s is %d',char(txt),hdr.nChans);

% begin sweep in ms
[txt,hdr.bsweep]=strread(char(hd{3}),'%s%f','delimiter','=');
fprintf('\n%s is %2.2f',char(txt),hdr.bsweep);

% sampling interval in ms
[txt,hdr.sampintms]=strread(char(hd{4}),'%s%f','delimiter','=');
fprintf('\n%s is %1.2f (or %2.1f Hz)',char(txt),hdr.sampintms,(1000./hdr.sampintms));
hdr.Fs = 1000./hdr.sampintms;

% bins per micro volt
[txt,hdr.binsuV]=strread(char(hd{5}),'%s%f','delimiter','=');
fprintf('\n%s is %1.2f',char(txt),hdr.binsuV);

% start time
tt=char(hd{6});
hdr.start_time=tt(end-7:end);
fprintf('\nStart Time is %s\n',hdr.start_time);

% header line 2: names of recording channels
hline2 = fgets(fid);

% channel names as cell array
%ch_names=textscan(hline2,'%s');
[hdr.label]=strread(hline2,'%s');
hdr.nTrials = 1;
hdr.nSamplesPre = 0;
% convert to char array
%ch_names=char(ch_names{1});


dc = textscan(fid,'%f',hdr.nChans*hdr.nSamples,'delimiter','\n','headerlines',0);
df = dc{1};clear dc;
df = df.*10; % multiply by 10, don't forget to divide when reading in data!
fprintf('converting ASCII floats to 16-bit signed ints');

dat = reshape(df,hdr.nChans,hdr.nSamples);
clear df;
dat = int16(dat)./10;
hdr.dat = double(dat);
clear dat;

% close input file
fclose(fid);