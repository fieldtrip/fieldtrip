function varargout = opensignals_txt(filename, hdr, begsample, endsample, chanindx)

% OPENSIGNALS_TXT reads time series data from a Bitalino OpenSignals txt file
%
% Use as
%   hdr = opensignals_txt(filename);
%   dat = opensignals_txt(filename, hdr, begsample, endsample, chanindx);
%   evt = opensignals_txt(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D

% Copyright (C) 2019 Robert Oostenveld
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

fid = fopen_or_error(fullname, 'r');
line1 = fgetl(fid);
line2 = fgetl(fid);
fclose(fid);
% convert the json string into a structure
orig = jsondecode(line2(2:end));
fn = fieldnames(orig);
if numel(fn)==1
  % it contains a single field with the date and time
  orig = orig.(fn{1});
end

dat = readtable(fullname, 'FileType', 'text');
keep = false(1,size(dat,2));
for i=1:size(dat,2)
  keep(i) = isnumeric(dat{:,i});
end
dat = table2array(dat(:,keep))';

if needhdr
  %% parse the header
  hdr.Fs = orig.samplingRate;
  hdr.nChans = length(orig.column);
  hdr.nTrials = 1;
  hdr.nSamplesPre = 0;
  hdr.nSamples = size(dat,2);
  hdr.label = orig.column;
  
  % remember the details
  hdr.orig = orig;
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% parse the data
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% parse the events
  ft_warning('cannot read events from %s', filename);
  
  varargout = {[]};
end
