function [trl, event] = ft_trialfun_brainvision_segmented(cfg)
% FT_TRIALFUN_BRAINVISION_SEGMENTED creates trials for a Brain Vision Analyzer
% dataset that was segmented in the BVA software.
%
% Use as 
%   cfg          = [];
%   cfg.dataset  = '<datasetname>.vhdr';
%   cfg.trialfun = 'ft_trialfun_brainvision_segmented';
%   cfg  = ft_definetrial(cfg);
%   data = ft_preprocessing(cfg);
%
% Optionally, you can also specify:
%   cfg.trigformat = 'S %d';
% which will instruct this function to parse stimulus triggers according to
% the format you specified. The default is 'S %d'. cfg.trigformat always
% needs to contain exactly one %d code. Trigger values read in this way
% will be stored in columns 4 and upwards of the output 'trl' matrix, and
% will end up in data.trialinfo if this matrix is subsequently passed to
% ft_preprocessing in order to read in data.
%
% A dataset needs to consist of three files: a .eeg, .vhdr, and .vmrk file.
% cfg.dataset needs to refer to the .vhdr file.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% Copyright (C) 20014, Robert Oostenveld, FCDC
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

hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% set the defaults
cfg.trigformat = ft_getopt(cfg, 'trigformat', 'S %d');

sel = find(strcmp({event.type}, 'New Segment'));
ntrial = length(sel);
fprintf('found %d segments in the data\n', ntrial);

sample = [event.sample];
if any(diff(sample)<0)
  % events should be ordered according to the sample in the file at which they occur
  warning('reordering events based on their sample number');
  [sample, order] = sort(sample);
  event = event(order);
end

% for Brain Vision Analyzer the event type and value are strings
type  = {event.type};
value = {event.value};

begsample = nan(1,ntrial);
endsample = nan(1,ntrial);
offset    = nan(1,ntrial);
stim      = cell(1,ntrial);

for i=1:ntrial
  begevt = sel(i);
  if i<ntrial
    endevt = sel(i+1);
  else
    endevt = length(event);
  end
  
  % specify at which data sample the trial begins
  begsample(i) = sample(sel(i));
  
  % find the sample on which the trial is time-locked
  t0 = find(strcmp(type((begevt+1):endevt), 'Time 0'));
  t0 = t0(1); % take the first one in case there are multiple
  offset(i) = -(sample(begevt+t0) - begsample(i));
  
  % construct a list of all stimulus events in each segment
  stim{i} = find(strcmp(type((begevt+1):endevt), 'Stimulus'))+begevt;
end

% the endsample of each trial aligns with the beginsample of the next one
endsample(1:end-1) = begsample(2:end);
% the last endsample corresponds to the end of the file
endsample(end)     = hdr.nSamples*hdr.nTrials;

% add the stimulus events to the output, if possible
numstim = cellfun(@length, stim);
if all(numstim==numstim(1))
  for i=1:length(stim)
    for j=1:numstim(i)
      stimvalue  = sscanf(value{stim{i}(j)}, cfg.trigformat);
      stimsample = sample(stim{i}(j));
      stimtime   = (stimsample - begsample(i) + offset(i))/hdr.Fs; % relative to 'Time 0'
      trialinfo(i,2*(j-1)+1) = stimvalue;
      trialinfo(i,2*(j-1)+2) = stimtime;
    end % j
  end % i
else
  warning('the trials have a varying number of stimuli, not adding them to the "trl" matrix');
  trialinfo = [];
end

% combine the sample information and the trial information
trl = cat(2, [begsample(:) endsample(:) offset(:)], trialinfo);

