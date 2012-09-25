function [hdr, dat] = read_neurosim_signals(filename)

% READ_NEUROSIM_SIGNALS reads the "signals" file that is written by Jan
% van der Eerden's NeuroSim software.
%
% See also FT_READ_HEADER, FT_READ_DATA

% Copyright (C) 2012 Robert Oostenveld
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

if isdir(filename)
  filename = fullfile(filename, 'signals');
end

needdat = (nargout>1);
needdat = true; % otherwise it fails to determine the number of samples in the file correctly

label = {};
orig  = {};

fid = fopen(filename, 'rb');

% read the header
line =  '#';
while ~isempty(line)
  % parse the content of the line, determine the label for each column
  colid = sscanf(line, '# column %d:', 1);
  if ~isempty(colid)
    label{colid} = strtrim(line(find(line==':')+1:end));
  end
  
  offset = ftell(fid); % remember the file pointer position
  line   = fgetl(fid); % get the next line
  if ~isempty(line) && line(1)~='#' && ~isempty(str2num(line))
    % the data starts here, rewind the last line
    fseek(fid, offset, 'bof');
    line = [];
  else
    orig{end+1} = line;
  end
end

if needdat
  % read the complete data
  dat = fscanf(fid, '%f', [length(label), inf]);
else
  % read only a small piece of data
  % this allows the code further down to determine the sampling frequency
  dat = fscanf(fid, '%f', [length(label), 10]);
end

fclose(fid);

% get the time axis, this is needed to determine the sampling frequency
time    = dat(match_str(label, 'Time'), :)/1e3; % this is in ms, not in seconds
fsample = median(1./diff(time));

% convert the header into fieldtrip style
hdr             = [];
hdr.label       = label(:);
hdr.Fs          = fsample;
hdr.nChans      = length(label);
% represent it as a single continuous recording
if needdat
  hdr.nSamples    = length(time);
else
  hdr.nSamples    = inf; % the total duration of the file is not known
end

hdr.FirstTimeStamp     = min(time);
hdr.TimeStampPerSample = 1e3/fsample; % how many timestamps (ms) are in one sample
hdr.nSamplesPre        = 0;
hdr.nTrials            = 1;

% also store the original ascii header details
hdr.orig        = orig(:);

