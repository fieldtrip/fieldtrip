function varargout = bucn_txt(filename, hdr, begsample, endsample, chanindx)

% BUCN_TXT reads the txt files produced by the UCL-Birkbeck NIRS machines, also known
% as the NTS fNIRS system from Gowerlabs. See https://www.gowerlabs.co.uk/nts
%
% See also READ_BUCN_NIRSHDR, READ_BUCN_NIRSDATA, READ_BUCN_NIRSEVENT, QUALISYS_TSV, MOTION_C3D

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

persistent data previous_fullname

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(previous_fullname) || ~isequal(fullname, previous_fullname) || isempty(data)
  % remember the full filename including path
  previous_fullname = fullname;

  % specifying false for ReadVariableNames fails on MATLAB R2015b, hence the 0
  % MATLAB R2015b does not have the 'MissingRule' option yet, you can specify 'HeaderLines' as 1 to skip the first line (or remove it manually from the txt file) 
  dat = readtable(filename, 'ReadVariableNames', 0, 'Delimiter', '\t', 'MissingRule', 'omitrow');
  dat = table2array(dat);

  [nrow, ncol] = size(dat);

  % each detector is a single row with (2 + nsource*2) values
  % each detector is sampled at a slightly different time

  ft_warning('assuming wavelengths of 780nm and 850nm (hard-coded)');
  wavelength = [780 850]; % FIXME

  nsource   = (size(dat, 2) - 2)/2; % minus time and trigger, divided by the number of wavelengths
  ndetector = sum((dat(1:100,1)-dat(1,1))<0.020); % look at the first 100 rows, consider them one sample if they are less than 20ms apart
  nsamples  = floor(nrow/ndetector);
  nchans    = nsource*ndetector*length(wavelength);

  %% construct the opto structure

  opto = [];
  opto.optolabel = {};
  opto.optotype = {};
  opto.optopos = nan(0,3);
  for i=1:nsource
    opto.optolabel{end+1} = sprintf('S%d', i);
    opto.optotype{end+1}  = 'transmitter';
    opto.optopos(end+1,:) = [nan nan nan];
  end
  for i=1:nsource
    opto.optolabel{end+1} = sprintf('D%d', i);
    opto.optotype{end+1}  = 'receiver';
    opto.optopos(end+1,:) = [nan nan nan];
  end
  opto.wavelength = wavelength;

  % construct all possible channels for all combinations of optodes and wavelengths
  selS = find(startsWith(opto.optolabel, 'S')); % sources or transmitters
  selD = find(startsWith(opto.optolabel, 'D')); % detectors or receivers
  labelS = opto.optolabel(selS);
  labelD = opto.optolabel(selD);

  nopto = length(selS) + length(selD);
  nchan = length(selS) * length(selD) * length(opto.wavelength);

  opto.label    = cell(nchan, 1);
  opto.chanpos  = nan(nchan, 3);
  opto.tra      = zeros(nchan, nopto);
  opto.unit     = 'mm';

  chan = 1;
  for i=1:length(labelS)
    for j=1:length(labelD)
      for k=1:length(wavelength)
        txIndex = selS(i); % transmitter/source
        rxIndex = selD(j); % receiver/detector
        opto.label{chan} = sprintf('%s-%s [%dnm]', labelS{i}, labelD{j}, opto.wavelength(k));

        txPos = opto.optopos(txIndex,:);
        rxPos = opto.optopos(rxIndex,:);
        opto.chanpos(chan,:) = (txPos + rxPos)/2; % place the channel in between the two optodes

        opto.tra(chan, txIndex) = +k;
        opto.tra(chan, rxIndex) = -k;
        chan = chan+1;
      end
    end
  end

  opto.label = opto.label(:);
  opto.optolabel = opto.optolabel(:);

  % there is only a single continuous segment
  data = [];
  data.time = {nan(1,nsamples)};
  data.trial = {nan(nchans+1,nsamples)}; % plus one for the trigger channel
  data.label = {};
  data.opto = opto;

  data.label{1} = 'rectime'; % add the recording time as first
  data.label{2} = 'trigger'; % add the trigger channel next
  for i=1:ndetector % for each detector, this is a single row
    for j=1:nsource % for each source, these are along columns, together with wavelengths
      data.label{end+1} = sprintf('S%d-D%d [%dnm]', j, i, wavelength(1));
      data.label{end+1} = sprintf('S%d-D%d [%dnm]', j, i, wavelength(2));
    end
  end
  data.label = data.label(:); % make it a column vector

  for sample=1:nsamples

    begrow = (sample-1)*ndetector + 1;
    endrow = begrow + ndetector - 1;

    tim = dat(begrow:endrow, 1); % time, in seconds
    trg = dat(begrow:endrow, 2); % trigger
    smp = dat(begrow:endrow, 3:end);

    rectim  = mean(tim); % take the mean as the time for all channels
    trigger = mode(trg); % take the most occuring value for the trigger

    % construct the time axis
    data.time{1}(sample)     = rectim;

    data.trial{1}(1, sample) = rectim;  % add the recording time as first
    data.trial{1}(2, sample) = trigger; % add the trigger channel next

    for i=1:ndetector % for each detector, this is a single row
      for j=1:(2*nsource) % for each source and wavelength, these are along columns
        chan = (i-1)*(2*nsource)+j;
        data.trial{1}(chan+2, sample) = smp(i,j);
      end
    end
  end % for samples

else
  % use the persistent data to speed up subsequent read operations
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if needhdr
  hdr = [];
  hdr.opto        = data.opto;
  hdr.label       = data.label(:);
  hdr.nChans      = length(hdr.label);
  hdr.Fs          = 1 / median(diff(data.time{1}));
  hdr.nSamples    = length(data.time{1});
  hdr.nSamplesPre = 0; % continuous data
  hdr.nTrials     = 1; % continuous data
  hdr.chantype    = repmat({'nirs'}, size(hdr.label));
  hdr.chantype{1} = 'trigger';
  hdr.chantype{2} = 'rectime';
  hdr.chanunit    = repmat({'unknown'}, size(hdr.label));

  % return the header
  varargout = {hdr};

elseif needdat
  dat = data.trial{1}(chanindx,begsample:endsample);

  % return the data
  varargout = {dat};

elseif needevt
  ft_warning('reading of events is not yet implemented');

  % return the events
  varargout = {[]};

end
