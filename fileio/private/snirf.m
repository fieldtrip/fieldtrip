function varargout = snirf(filename, hdr, begsample, endsample, chanindx)

% SNIRF reads data from a SNIRF file and returns it in a format that FieldTrip understands.
% See https://github.com/fNIRS
%
% Use as
%   hdr = snirf(filename);
%   dat = snirf(filename, hdr, begsample, endsample, chanindx);
%   evt = snirf(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D

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

% this uses the SNIRF reading functions from the Homer3 toolbox
ft_hastoolbox('homer3', 1);

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

snirf = SnirfClass;
snirf.Load(filename);

hasdata = ~isempty(snirf.data) && ~isempty(snirf.data.dataTimeSeries);
hasaux = ~isempty(snirf.aux)   && ~isempty(snirf.aux.dataTimeSeries);

if      hasdata &&  hasaux
  ndata    = length(snirf.data.measurementList);
  naux     = min(size(snirf.aux.dataTimeSeries)); % assume there are more samples than AUX channels
  nsamples = length(snirf.data.time);
  assert(isequal(length(snirf.aux.time), nsamples), 'the number of samples does not match');
elseif ~hasdata &&  hasaux
  ndata    = 0;
  naux     = min(size(snirf.aux.dataTimeSeries)); % assume there are more samples than AUX channels
  nsamples = length(snirf.aux.time);
elseif  hasdata && ~hasaux
  ndata    = length(snirf.data.measurementList);
  naux     = 0;
  nsamples = length(snirf.data.time);
else
  ft_error('this file contains no data');
end

if hasdata
  % the code below expects the timeseries data to be organized as nsamples*nchans
  if ~isequal(size(snirf.data.dataTimeSeries), [nsamples, ndata])
    snirf.data.dataTimeSeries = snirf.data.dataTimeSeries';
  end
end

if hasaux
  % the code below expects the timeseries data to be organized as nsamples*nchans
  if ~isequal(size(snirf.aux.dataTimeSeries), [nsamples, naux])
    snirf.aux.dataTimeSeries = snirf.aux.dataTimeSeries';
  end
end

if needhdr
  %% read the header
  
  hdr.Fs       = 1/median(diff(snirf.data.time));
  hdr.nChans   = ndata + naux;
  hdr.nSamples = nsamples;
  hdr.nSamplesPre = 0;  % assume it is continuous
  hdr.nTrials  = 1;     % assume it is continuous
  hdr.label    = {};
  hdr.chantype = {};
  hdr.chanunit = {};
  for i=1:ndata
    % construct the channel name using the source, detector and wavelength
    try
      d = snirf.probe.detectorLabels{snirf.data.measurementList(i).detectorIndex}; % receiver
      s = snirf.probe.sourceLabels{snirf.data.measurementList(i).sourceIndex};     % transmitter
      w = snirf.probe.wavelengths(snirf.data.measurementList(i).wavelengthIndex);
      hdr.label   {end+1} = sprintf('%s-%s [%dnm]', d, s, w);
    catch
      % it is apparently possible that the probe information is not specified
      d = snirf.data.measurementList(i).detectorIndex;   % receiver
      s = snirf.data.measurementList(i).sourceIndex;     % transmitter
      w = snirf.data.measurementList(i).wavelengthIndex;
      hdr.label   {end+1} = sprintf('D%d-S%d [%d]', d, s, w);
    end
    hdr.chantype{end+1} = 'nirs';
    hdr.chanunit{end+1} = 'unknown';
  end
  for i=1:naux
    hdr.label   {end+1} = sprintf('aux%d', i);
    hdr.chantype{end+1} = 'aux';
    hdr.chanunit{end+1} = 'unknown';
  end
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% read the data
  
  % concatenate the different types of data
  if      hasdata &&  hasaux
    dat = cat(1, snirf.data.dataTimeSeries', snirf.aux.dataTimeSeries');
  elseif ~hasdata &&  hasaux
    dat = snirf.aux.dataTimeSeries';
  elseif  hasdata && ~hasaux
    dat = snirf.data.dataTimeSeries';
  end
  
  % make selection of channels and samples
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% read the events
  
  ft_warning('reading events from a SNIRF file is not yet implemented');
  evt = [];
  
  % return the events
  varargout = {evt};
  
end
