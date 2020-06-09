function varargout = openbci_txt(filename, hdr, begsample, endsample, chanindx)

% OPENBCI_TXT reads time series data from a OpenBCI txt file
%
% Use as
%   hdr = openbci_txt(filename);
%   dat = openbci_txt(filename, hdr, begsample, endsample, chanindx);
%   evt = openbci_txt(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D, OPENSIGNALS_TXT

% Copyright (C) 2020 Robert Oostenveld
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

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

% use the full filename including path to distinguish between similarly named files in different directories
[p, f, x] = fileparts(filename);
if isempty(p)
  % no path was specified
  fullname = which(filename);
elseif startsWith(p, ['.' filesep])
  % a relative path was specified
  fullname = fullfile(pwd, p(3:end), [f, x]);
else
  fullname = filename;
end

orig = [];

fid = fopen_or_error(fullname, 'r');
line = '%';
while ~feof(fid) && line(1)=='%'
  line = fgetl(fid);
  if contains(line, '=')
    tok = strsplit(line, '=');
    key = lower(strip(tok{1}(2:end)));
    key(key==' ') = '_';
    val = strip(tok{2});
    switch key
      case 'number_of_channels'
        orig.(key) = str2double(val);
      case 'sample_rate'
        orig.(key) = str2double(strtok(val));
      otherwise
        orig.(key) = val;
    end
  end
end

fclose(fid);

% read the data as table
tab = readtable(fullname, 'FileType', 'text');

% remove all non-numeric columns
sel = cellfun(@isnumeric, table2cell(tab(1,:)));
% also remove the 1st column with the SampleIndex
sel(1) = false;
% select columns of interest and transpose
dat = table2array(tab(:,sel))';

if needhdr
  %% construct the header
  hdr.Fs = orig.sample_rate;
  hdr.nChans = size(dat,1);
  hdr.nTrials = 1;
  hdr.nSamplesPre = 0;
  hdr.nSamples = size(dat,2);
  hdr.label = {};
  hdr.chantype = {};
  hdr.chanunit = {};
  for i=1:hdr.nChans
    hdr.label{i} = num2str(i);
    hdr.chantype{i} = 'unknown';
    hdr.chanunit{i} = 'unknown';
  end
  hdr.chantype(1:orig.number_of_channels) = {'eeg'};
  
  % remember the details
  hdr.orig = orig;
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% select the requested data
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% parse the events
  ft_warning('cannot read events from %s', filename);
  
  varargout = {[]};
end
