function [hdr] = read_bucn_nirshdr(filename)

% READ_BUCN_NIRSHDR reads the header information of ASCII-formatted NIRS
% data acquired with the UCL-BIRKBECK machine and postprocessed by the
% Paris group. The first line contains the channel labels and the rest of
% the file contains per line a time sample. The first column specifies the
% time axis.
%
% Use as
%   [hdr] = read_bucn_nirshdr(filename)
%
% See also read_bucn_nirsdata, READ_BUCN_NIRSEVENT

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
% $Id: read_bucn_nirshdr.m$

fid = fopen(filename, 'r');

% read the first line
line1 = textscan(fid, '%[^\n]',1);

% field delimiter can be space or tab
labelspc = textscan(line1{1}{1}, '%[^ ]');
labeltab = textscan(line1{1}{1}, '%[^\t]');

% let tab as a delimiter prevail
if numel(labeltab{1})>1
  label = labeltab{1};
else
  label = labelspc{1};
end
nchan = numel(label);
Fs    = str2num(strtok(strtok(label{1},'#Time.'),'Hz'));

% test whether the channel labels are non-numeric
labelnumber = cellfun(@str2num, label, 'UniformOutput', false);
labelstring = cellfun(@isempty, labelnumber, 'UniformOutput', true);
if ~any(labelstring)
  ft_error('channel labels were not found in the first line of the file');
end

% read the rest
dat = textscan(fid, '%f');
fclose(fid);

dat  = reshape(dat{1}, nchan, []);
nsmp = size(dat,2);

% create the output
hdr          = [];
hdr.Fs       = Fs;
hdr.label    = label;
hdr.nTrials  = 1;
hdr.nSamples = nsmp;
hdr.nSamplesPre = 0;
hdr.nChans   = nchan;
hdr.time     = dat(1,:); % events in the raw event file have both a sample and a time stamp
