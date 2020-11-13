function varargout = snirf(filename, hdr, begsample, endsample, chanindx)

% SNIRF reads data from a SNIRF file and returns it in a format that FieldTrip understands.
%
% See https://github.com/fNIRS/snirf/blob/master/snirf_specification.md
%
% Use as
%   hdr = snirf(filename);
%   dat = snirf(filename, hdr, begsample, endsample, chanindx);
%   evt = snirf(filename, hdr);
%
% The SNIRF format allows for multiple blocks of data channels anx aux channels, each
% with a different sampling frequency. That is not allowed in this code; all channels
% must have the same sampling rate and be sampled at the same time.
%
% See also SNIRF2OPTO, FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D

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

% this uses the SNIRF reading functions from the Homer3 toolbox
ft_hastoolbox('homer3', 1);

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

snirf = SnirfClass; % this is from Homer3
snirf.Load(filename);

hasdata = ~isempty(snirf.data) && ~isempty(snirf.data(1).dataTimeSeries);
hasaux = ~isempty(snirf.aux)   && ~isempty(snirf.aux(1).dataTimeSeries);

ndata = 0;
if hasdata
  for i=1:length(snirf.data)
    ndata = ndata + length(snirf.data(i).measurementList);
    assert(isequal(snirf.data(i).time, snirf.data(1).time), 'different sampling rates are not supported');
  end
else
  % code it as empty, not as 0x1 or 1x0
  snirf.data = [];
end

naux = 0;
if hasaux
  for i=1:length(snirf.aux)
    naux = naux + 1; % each aux channel is stored separately
    assert(isequal(snirf.aux(i).time, snirf.aux(1).time), 'different sampling rates are not supported');
  end
else
  % code it as empty, not as 0x1 or 1x0
  snirf.aux = [];
end

if hasdata && hasaux
  assert(isequal(snirf.data(1).time, snirf.aux(1).time), 'different sampling rates are not supported');
  nsamples = length(snirf.data(1).time);
elseif hasdata
  nsamples = length(snirf.data(1).time);
elseif hasaux
  nsamples = length(snirf.aux(1).time);
else
  ft_error('incorrect data format');
end

if hasdata
  % the code below expects the timeseries data to be organized as nsamples*nchans
  for i=1:length(snirf.data)
    if ~isequal(size(snirf.data(i).dataTimeSeries), [nsamples, ndata])
      snirf.data(i).dataTimeSeries = snirf.data(i).dataTimeSeries';
    end
  end
end

if hasaux
  % the code below expects the timeseries data to be organized as nsamples*nchans
  for i=1:length(snirf.aux)
    if ~isequal(size(snirf.aux(i).dataTimeSeries), [nsamples, 1])
      snirf.aux(i).dataTimeSeries = snirf.aux(i).dataTimeSeries';
    end
  end
end

if needhdr
  %% read the header
  
  if hasdata
    hdr.Fs     = 1/median(diff(snirf.data(1).time));
  elseif hasaux
    hdr.Fs     = 1/median(diff(snirf.aux(1).time));
  end
  hdr.nChans   = ndata + naux;
  hdr.nSamples = nsamples;
  hdr.nSamplesPre = 0;  % assume it is continuous
  hdr.nTrials  = 1;     % assume it is continuous
  hdr.label    = {};
  hdr.chantype = {};
  hdr.chanunit = {};
  for i=1:length(snirf.data)
    for j=1:length(snirf.data(i).measurementList)
      % construct the channel name using the source, detector and wavelength
      try
        d = snirf.probe.detectorLabels{snirf.data(i).measurementList(j).detectorIndex};     % receiver
        s = snirf.probe.sourceLabels  {snirf.data(i).measurementList(j).sourceIndex};       % transmitter
        w = snirf.probe.wavelengths   (snirf.data(i).measurementList(j).wavelengthIndex);
        hdr.label   {end+1} = sprintf('%s-%s [%dnm]', s, d, round(w));
      catch
        % it is apparently possible that the probe information is not specified
        d = snirf.data(i).measurementList(j).detectorIndex;   % receiver
        s = snirf.data(i).measurementList(j).sourceIndex;     % transmitter
        w = snirf.data(i).measurementList(j).wavelengthIndex;
        hdr.label   {end+1} = sprintf('S%d-D%d [%d]', s, d, w);
      end
      hdr.chantype{end+1} = 'nirs';
      hdr.chanunit{end+1} = 'unknown';
    end
  end
  
  for i=1:length(snirf.aux)
    hdr.label   {end+1} = snirf.aux(i).name;
    hdr.chantype{end+1} = 'aux';
    hdr.chanunit{end+1} = 'unknown';
  end
  
  % convert the probe and measurementList to a FieldTrip opto structure
  try
    hdr.opto = snirf2opto(snirf.probe, snirf.data.measurementList);
  catch
    ft_warning('SNIRF probe and measurementList are inconsistent');
  end
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% read the data
  % https://github.com/fNIRS/snirf/blob/master/snirf_specification.md#nirsidataj
  
  % concatenate the different types of data
  if      hasdata &&  hasaux
    dat = {snirf.data.dataTimeSeries};
    aux = {snirf.aux.dataTimeSeries};
  elseif ~hasdata &&  hasaux
    dat = {};
    aux = {snirf.aux.dataTimeSeries};
  elseif  hasdata && ~hasaux
    dat = {snirf.data.dataTimeSeries};
    aux = {};
  end
  dat = cat(2, dat{:}, aux{:})';
  
  % make selection of channels and samples
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% read the events
  % see https://github.com/fNIRS/snirf/blob/master/snirf_specification.md#nirsistimj
  
  evt = [];
  
  try
    for i=1:numel(snirf.stim)
      for j=1:size(snirf.stim(i).data,1)
        evt(end+1).type      = snirf.stim(i).name;
        evt(end  ).sample    = round(snirf.stim(i).data(j,1) * hdr.Fs + 1); % time 0 corresponds to sample 1
        evt(end  ).duration  = round(snirf.stim(i).data(j,2) * hdr.Fs);
        evt(end  ).value     = snirf.stim(i).data(j,3);
      end
    end
    % sort the events on their sample number, this is consistent with FT_READ_EVENT
    [~, indx] = sort([evt.sample]);
    evt = evt(indx);
  catch
    % there are no events
  end
  
  % return the events
  varargout = {evt};
  
end
