function [trl] = ft_trialfun_twoclass_classification(cfg)

% FT_TRIALFUN_TWOCLASS_CLASSIFICATION can be used to train and test a real-time
% classifier in offline and online mode. It selects pieces of data in the two classes
% based on two trigger values. The first N occurences in each class are marked as
% training items. All subsequent occurrences are marked as test items.
%
% This function can be used in conjunction with FT_REALTIME_CLASSIFICATION. The
% configuration structure should contain
%   cfg.dataset              = string with the filename
%   cfg.trialfun             = 'ft_trialfun_twoclass_classification'
%   cfg.trialdef.numtrain    = number of training items, e.g. 20
%   cfg.trialdef.eventvalue1 = trigger value for the 1st class
%   cfg.trialdef.eventvalue2 = trigger value for the 2nd class
%   cfg.trialdef.eventtype   = string, e.g. 'trigger'
%   cfg.trialdef.prestim     = latency in seconds, e.g. 0.3
%   cfg.trialdef.poststim    = latency in seconds, e.g. 0.7
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL

% Copyright (C) 2009, DCCN
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

% these are used to count the number of training items in each class
persistent numtrain1
persistent numtrain2

if isempty(numtrain1)
  numtrain1 = 0;
end
if isempty(numtrain2)
  numtrain2 = 0;
end

if isfield(cfg, 'hdr')
  hdr = cfg.hdr;
else
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
end

if isfield(cfg, 'event')
  event = cfg.event;
else
  event = ft_read_event(cfg.headerfile, 'headerformat', cfg.headerformat);
end

baseline = round(cfg.trialdef.prestim*hdr.Fs);
duration = round(cfg.trialdef.prestim*hdr.Fs + cfg.trialdef.poststim*hdr.Fs);

% make a subset of the interesting events
sel   = strcmp(cfg.trialdef.eventtype, {event.type});
event = event(sel);
num   = length(event);
trl   = zeros(num,5);

for i=1:num
  % determine the location of this trial in the data stream
  begsample = event(i).sample - baseline;
  endsample = begsample + duration - 1;
  offset    = baseline;
  % determine the class and wether this trial is eligeable for training
  if event(i).value==cfg.trialdef.eventvalue1
    class = 1;
    train = (numtrain1 < cfg.trialdef.numtrain); % boolean
    numtrain1 = numtrain1 + train;  % increment the counter
  elseif event(i).value==cfg.trialdef.eventvalue2
    class = 2;
    train = (numtrain2 < cfg.trialdef.numtrain); % boolean
    numtrain2 = numtrain2 + train;  % increment the counter
  else
    % the class is unknown and therefore irrelevant
    class = nan;
    train = false;
  end
  % remember this trial, the class and whether it should be used for training
  trl(i,:) = [begsample endsample offset class train];
end
