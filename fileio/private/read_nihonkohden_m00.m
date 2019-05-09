function [hdr, dat] = read_nihonkohden_m00(filename)

% READ_NIHONKOHDEN_M00 reads the header and data from a file in the Nihon Kohden *.m00 format.
% This implementation is an adaptation of convert_nkascii2mat.m and get_nkheader.m written
% by Timothy Ellmore, see https://openwetware.org/wiki/Beauchamp:AnalyzeEEGinMatlab.
%
% Use as
%   [hdr, dat] = read_nihonkohden_m00(filename)
%
% This returns a FieldTrip compatible header structure and the data matrix.
%
% See also FT_READ_HEADER, FT_READ_DATA

% Copyright (C) 2007, Timothy.Ellmore@uth.tmc.edu
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

% open the ascii file
fid = fopen_or_error(filename, 'r');

%header line 1: acquisition information
hline1 = fgets(fid);

% get the six different fields in the header, store in cell-array
hd = textscan(hline1, '%s %s %s %s %s %s');
fprintf('\nNihon Kohden ASCII EEG header fields:');
fprintf('\n-------------------------------------');

% number of timepoints
[txt, ntpoints] = strread(char(hd{1}), '%s%d', 'delimiter', '=');
fprintf('\n%s is %d', char(txt), ntpoints);

% number of channels sampled
[txt, nchannels] = strread(char(hd{2}), '%s%d', 'delimiter', '=');
fprintf('\nNumber of %s is %d', char(txt), nchannels);

% begin sweep in ms
[txt, bsweep] = strread(char(hd{3}), '%s%f', 'delimiter', '=');
fprintf('\n%s is %2.2f', char(txt), bsweep);

% sampling interval in ms
[txt, sampintms] = strread(char(hd{4}), '%s%f', 'delimiter', '=');
fprintf('\n%s is %1.2f (or %2.1f Hz)', char(txt), sampintms, (1000./sampintms));

% bins per micro volt
[txt, binsuV] = strread(char(hd{5}), '%s%f', 'delimiter', '=');
fprintf('\n%s is %1.2f', char(txt), binsuV);

% start time
tt = char(hd{6});
start_time = tt(end-7:end);
fprintf('\nStart Time is %s\n', start_time);

% header line 2: names of recording channels
hline2 = fgets(fid);

% channel names as cell-array
ch_names = textscan(hline2, '%s');

% convert it into a FieldTrip header
hdr.Fs          = 1000./sampintms;
hdr.nChans      = nchannels;
hdr.nSamples    = ntpoints;
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;
hdr.label       = ch_names{1};

% keep the original header details
hdr.orig.hline1     = hline1;
hdr.orig.hline2     = hline2;
hdr.orig.binsuV     = binsuV;
hdr.orig.bsweep     = bsweep;
hdr.orig.ntpoints   = ntpoints;
hdr.orig.sampintms  = sampintms;
hdr.orig.start_time = start_time;

if nargout>1
  % read the remaining data
  dat = textscan(fid, '%f', hdr.nChans*hdr.nSamples);
  dat = reshape(dat{1}, hdr.nChans, hdr.nSamples);
  % calibrate the data to obtain uV values
  dat = dat ./ hdr.orig.binsuV;
end

% close input file
fclose(fid);
