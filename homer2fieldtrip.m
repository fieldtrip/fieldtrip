function [data, event] = homer2fieldtrip(nirs, varargin)

% HOMER2FIELDTRIP converts a continuous raw data structure from Homer to FieldTrip
% format.
%
% Use as
%   data = homer2fieldtrip(filename)
% where the input is a file name, or
%   data = homer2fieldtrip(nirs)
% where the input nirs structure is according to the Homer format and the output data
% structure is formatted according to the output of FT_PREPROCESSING.
%
% See https://www.nitrc.org/plugins/mwiki/index.php/homer2:Homer_Input_Files#NIRS_data_file_format
% for a description of the Homer data structure.
%
% See also FIELDTRIP2HOMER, FT_PREPROCESSING, FT_DATATYPE_RAW

% Copyright (C) 2020, Robert Oostenveld
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

% this allows to select channels and to override the automatic trigger channel detection
chanindx = ft_getopt(varargin, 'chanindx');

if ischar(nirs)
  % Homer files are MATLAB files in disguise
  % see https://www.nitrc.org/plugins/mwiki/index.php/homer2:Homer_Input_Files#NIRS_data_file_format
  nirs = load(nirs, '-mat');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the data similar to FT_READ_DATA

% copy the data from the raw intensity time courses
dat = nirs.d;
if isfield(nirs, 'aux')
  % concatenate the AUX channel at the end, consistent with FT_READ_HEADER
  dat = [dat nirs.aux];
end
if isfield(nirs, 's')
  % concatenate the stimulus channel at the end, consistent with FT_READ_HEADER
  dat = [dat nirs.s];
end

% transpose the data, FieldTrip expects nchan*ntime
dat = dat';

% the Homer structure includes the time axis, ensure it is a row vector
time = nirs.t(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the header similar to FT_READ_HEADER

hdr.label       = {};
hdr.nChans      = size(nirs.d,2);
hdr.nSamples    = size(nirs.d,1);
hdr.nSamplesPre = 0;
hdr.nTrials     = 1; % assume continuous data, not epoched
hdr.Fs          = 1/median(diff(nirs.t));

% number of wavelengths times sources times detectors
if ~isfield(nirs.SD, 'nSrcs')
  nirs.SD.nSrcs = size(nirs.SD.SrcPos,1);
end
if ~isfield(nirs.SD, 'nDets')
  nirs.SD.nDets = size(nirs.SD.DetPos,1);
end
assert(numel(nirs.SD.Lambda)*nirs.SD.nSrcs*nirs.SD.nDets >= hdr.nChans);

try
  % use the transmitter and receiver numbers and the wavelength to form the the channel names
  for i=1:hdr.nChans
    tx = nirs.SD.MeasList(i,1); % transmitter
    rx = nirs.SD.MeasList(i,2); % receiver
    wl = nirs.SD.Lambda(nirs.SD.MeasList(i,4)); % wavelength in nm
    hdr.label{i} = sprintf('S%d-D%d [%dnm]', tx, rx, round(wl));
  end
catch
  ft_warning('creating default channel names');
  for i=1:hdr.nChans
    hdr.label{i} = num2str(i);
  end
end

hdr.chantype = repmat({'nirs'}, hdr.nChans, 1);
hdr.chanunit = repmat({'unknown'}, hdr.nChans, 1);

if isfield(nirs, 'aux')
  % concatenate the AUX channel(s) at the end, consistent with FT_READ_DATA
  if size(nirs.aux,2)==1
    hdr.nChans = hdr.nChans + 1;
    hdr.label{end+1} = 'aux';
    hdr.chantype{end+1} = 'aux';
    hdr.chanunit{end+1} = 'unknown';
  else
    for i=1:size(nirs.aux,2)
      hdr.nChans = hdr.nChans + 1;
      hdr.label{end+1} = sprintf('aux%d', i);
      hdr.chantype{end+1} = 'aux';
      hdr.chanunit{end+1} = 'unknown';
    end
  end % if single or multiple
end % if aux channels present

if isfield(nirs, 's')
  
  % the CondNames field is not always present, it is documented on page 28 in http://www.nmr.mgh.harvard.edu/martinos/software/homer/HOMER2_UsersGuide_121129.pdf
  if ~isfield(nirs, 'CondNames')
    if size(nirs.s,2)==1
      % do not include the column number, there is only one
      nirs.CondNames = {'s'};
    else
      for i=1:size(nirs.s,2)
        % include the column number
        nirs.CondNames{i} = sprintf('s%d', i);
      end
    end % if single or multiple
  end
  
  % concatenate the stimulus channel(s) at the end, consistent with FT_READ_DATA
  for i=1:size(nirs.s,2)
    hdr.nChans = hdr.nChans + 1;
    hdr.label{end+1} = nirs.CondNames{i};
    hdr.chantype{end+1} = 'stimulus';
    hdr.chanunit{end+1} = 'unknown';
  end
  
end % if stimulus channels present

% convert the measurement configuration details to an optode structure
hdr.opto = homer2opto(nirs.SD);

% keep all details except the data
hdr.orig = removefields(nirs, {'d', 't', 's', 'aux'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the FieldTrip raw data structure, see also FT_DATATYPE_RAW

data = [];
data.opto = hdr.opto;
data.hdr = hdr;
data.fsample = hdr.Fs;
data.sampleinfo = [1 hdr.nSamples*hdr.nTrials];
% the output has a single trial, representing continuous data
if isempty(chanindx)
  data.label = hdr.label;
  data.trial{1} = dat;
  data.time{1} = time;
else
  data.label = hdr.label(chanindx);
  data.trial{1} = dat(chanindx,:);
  data.time{1} = time;
end

% this should be a column vector
data.label = data.label(:);

% ensure that the output is according to the most recent standards
data = ft_checkdata(data, 'datatype', 'raw');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct an event structure similar to FT_READ_EVENT

if nargout>1 && isfield(nirs, 's')
  event = boolvec2event(nirs.s', 'type', nirs.CondNames);
end
