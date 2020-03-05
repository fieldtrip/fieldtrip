function [trl, event] = ft_trialfun_brainvision_segmented(cfg)

% FT_TRIALFUN_BRAINVISION_SEGMENTED creates trials for a Brain Vision Analyzer
% dataset that was segmented in the BVA software.
%
% Use as 
%   cfg          = [];
%   cfg.dataset  = 'filename.vhdr';
%   cfg.trialfun = 'ft_trialfun_brainvision_segmented';
%   cfg  = ft_definetrial(cfg);
%   data = ft_preprocessing(cfg);
%
% Optionally, you can specify:
%   cfg.stimformat = 'S %d'
% 
% The stimformat instruct this function to parse stimulus triggers according to
% the specific format. The default is 'S %d'. The cfg.stimformat always
% needs to contain exactly one %d code. The trigger values parsed in this way
% will be stored in columns 4 and upwards of the output 'trl' matrix, and
% after FT_PREPROCESSING will end up in data.trialinfo.
%
% A BrainVision dataset consists of three files: an .eeg, .vhdr, and a .vmrk 
% file. The option cfg.dataset should refer to the .vhdr file.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% Copyright (C) 2014, Robert Oostenveld, DCCN
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

% for backward compatibility
cfg = ft_checkconfig(cfg, 'renamed', {'trigformat', 'stimformat'});

% set the defaults
cfg.stimformat = ft_getopt(cfg, 'stimformat', 'S %d');

sel = find(strcmp({event.type}, 'New Segment'));
ntrial = length(sel);
fprintf('found %d segments in the data\n', ntrial);

sample = [event.sample];
if any(diff(sample)<0)
  % events should be ordered according to the sample in the file at which they occur
  ft_warning('reordering events based on their sample number');
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
endsample(1:end-1) = begsample(2:end) - 1;
% the last endsample corresponds to the end of the file
endsample(end)     = hdr.nSamples*hdr.nTrials;

% add the stimulus events to the output, if possible
numstim = cellfun(@length, stim);
if (numstim > 0)
  if all(numstim==numstim(1))
    for i=1:length(stim)
      for j=1:numstim(i)
        if isempty(value{stim{i}(j)})
          ft_error('missing a stimulus type definition in the related *.vmrk file');
        end
        stimvalue  = sscanf(value{stim{i}(j)}, cfg.stimformat);
        stimsample = sample(stim{i}(j));
        stimtime   = (stimsample - begsample(i) + offset(i))/hdr.Fs; % relative to 'Time 0'
        if isempty(stimvalue)
          ft_error('the upper case letter of the stimulus value does not match with definition of "cfg.stimformat"'); 
        end
        trialinfo(i,2*(j-1)+1) = stimvalue;
        trialinfo(i,2*(j-1)+2) = stimtime;
      end % j
    end % i
  else
    ft_warning('the trials have a varying number of stimuli, not adding them to the "trl" matrix');
    trialinfo = [];
  end
else
  ft_warning('the trials have no stimuli, no "trialinfo" will be added to the "trl" matrix');
  trialinfo = [];
end

% combine the sample information and the trial information
trl = cat(2, [begsample(:) endsample(:) offset(:)], trialinfo);

