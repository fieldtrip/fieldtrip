function [trl, allevents] = ft_trialfun_bids(cfg)

% FT_TRIALFUN_BIDS determines trials/segments to be used for subsequent analysis, on
% the basis of the BIDS "events.tsv" file. This function should in general not be
% called directly, it will be called by FT_DEFINETRIAL.
%
% Use this function by calling 
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset   = string with the filename
%   cfg.trialdef  = structure with the details of trial definition, see below
%   cfg.trialfun  = 'ft_trialfun_bids'
%
% The trialdef structure should contain the following specifications
%   cfg.trialdef.prestim    = latency in seconds (required)
%   cfg.trialdef.poststim   = latency in seconds (required)
% and you can specify your selection of events as
%   cfg.trialdef.columnname = columnvalue
% where the column name and value have to match those present in the events.tsv file.
%
% For example
%   cfg.trialdef.prestim  = 0.2;
%   cfg.trialdef.poststim = 0.8;
%   cfg.trialdef.task     = 'notarget';
%   cfg.trialdef.category = 'tools';
%   cfg.trialdef.modality = {'written', 'spoken'};
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL

% Copyright (C) 2021, Robert Oostenveld
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

% most defaults are in trialdef
cfg.trialdef = ft_getopt(cfg, 'trialdef', struct());

% specify the default file formats
cfg.eventformat   = ft_getopt(cfg, 'eventformat');
cfg.headerformat  = ft_getopt(cfg, 'headerformat');
cfg.dataformat    = ft_getopt(cfg, 'dataformat');

% get the header, this is among others for the sampling frequency
if isfield(cfg, 'hdr')
  ft_info('using the header from the configuration structure\n');
  hdr = cfg.hdr;
else
  % read the header, contains the sampling frequency
  ft_info('reading the header from ''%s''\n', cfg.headerfile);
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
end

% get the events
if isfield(cfg, 'event')
  ft_info('using the events from the configuration structure\n');
  events = cfg.event;
else
  % do not use FT_READ_EVENTS, as that will force the events in a structure
  eventsfile = bids_sidecar(cfg.dataset, 'events');
  events = ft_read_tsv(eventsfile);
end

% these cannot be used as column headings, as they are used as normal cfg options
assert(~ismember('prestim', events.Properties.VariableNames));
assert(~ismember('poststim', events.Properties.VariableNames));

% make a selection of the rows
sel = true(size(events,1),1);

fn = fieldnames(cfg.trialdef);
for i=1:numel(fn)
  if ismember(fn{i}, events.Properties.VariableNames)
    columnname  = fn{i};
    columnvalue = cfg.trialdef.(columnname);
    sel = sel & ismember(events.(columnname), columnvalue);
  end
end

% remember and return all events
allevents = events;
% continue with the selection of events
events = events(sel,:);

if ismember('sample', events.Properties.VariableNames)
  begsample = round(events.sample - cfg.trialdef.prestim * hdr.Fs);
  endsample = round(events.sample + cfg.trialdef.poststim * hdr.Fs);
  offset    = round(-cfg.trialdef.prestim*hdr.Fs) * ones(size(begsample));
else
  % use the onset, which is in seconds
  begsample = round((events.onset * hdr.Fs) + 1 - cfg.trialdef.prestim);
  endsample = round((events.onset * hdr.Fs) + 1 + cfg.trialdef.poststim);
  offset    = round(-cfg.trialdef.prestim*hdr.Fs) * ones(size(begsample));
end

trl = table(begsample, endsample, offset);
trl = cat(2, trl, events); % keep the selected event details

% remove trials that fall outside the data
outside = trl.begsample<1 | trl.endsample>(hdr.nSamples*hdr.nTrials);
if any(outside)
  ft_notice('removing %d trials that extend outside the data', sum(outside));
  trl = trl(~outside, :);
end
