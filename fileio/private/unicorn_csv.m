function varargout = unicorn_csv(filename, hdr, begsample, endsample, chanindx)

% UNICORN_CSV reads EEG data from the Gtec/Unicorn Hybrid Black
%
% See http://unicorn-bi.com/
%
% Use as
%   hdr = unicorn_csv(filename);
%   dat = unicorn_csv(filename, hdr, begsample, endsample, chanindx);
%   evt = unicorn_csv(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, MOTION_C3D, QUALISYS_TSV, LIBERTY_CSV

% Copyright (C) 2022, Robert Oostenveld
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

persistent csv previous_fullname

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

if isempty(previous_fullname) || ~isequal(fullname, previous_fullname) || isempty(csv)
  % remember the full filename including path
  previous_fullname = fullname;
  % read and remember the file content
  ws = warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
  csv = readtable(fullname, 'ReadVariableNames', true);
  warning(ws); % revert to the previous warning state
else
  % use the persistent variable to speed up subsequent read operations
end

if needhdr
  %% parse the header
  hdr = [];
  hdr.Fs          = 250; % the sampling rate is fixed
  hdr.nChans      = size(csv,2);
  hdr.nSamples    = size(csv,1);
  hdr.nSamplesPre = 1;
  hdr.nTrials     = 1;
  hdr.label       = csv.Properties.VariableNames';
  hdr.chantype    = repmat({'unknown'}, size(hdr.label));
  hdr.chanunit    = repmat({'unknown'}, size(hdr.label));

  hdr.chantype(startsWith(hdr.label, 'EEG'))            = {'eeg'};
  hdr.chantype(startsWith(hdr.label, 'Accelerometer'))  = {'accel'};
  hdr.chantype(startsWith(hdr.label, 'Gyroscope'))      = {'gyro'};

  % return the header details
  varargout = {hdr};

elseif needdat
  % make the selection of channels and samples
  dat = table2array(csv(begsample:endsample, chanindx))';

  % return the data
  varargout = {dat};

elseif needevt
  %% there are no events
  event = [];

  % return the events
  varargout = {event};

end
