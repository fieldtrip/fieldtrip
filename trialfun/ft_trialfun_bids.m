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
% The trialdef structure should either contain the following
%   cfg.trialdef.prestim    = latency in seconds
%   cfg.trialdef.poststim   = latency in seconds
% or the duration and offset relative to the event of interest
%   cfg.trialdef.duration    = latency in seconds
%   cfg.trialdef.offset      = latency in seconds
%
% You can specify your selection of events as
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

% Copyright (C) 2021-2024, Robert Oostenveld
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
  % do not use FT_READ_EVENT, as that will force the events in a structure
  eventsfile = bids_sidecar(cfg.dataset, 'events');
  events = ft_read_tsv(eventsfile);
end

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

% continue defining trials with the selection of events
events = events(sel,:);

% FIXME it would also be possible to construct trials if the events contained prestim and poststim
% FIXME it would also be possible to construct trials if the events contained duration and offset
% FIXME in those cases the user-specified cfg.trialdef options could overrule the ones from the events

if isfield(cfg.trialdef, 'prestim') && isfield(cfg.trialdef, 'poststim')
  % these are mutually exclusive
  assert(~isfield(cfg.trialdef, 'duration'));
  assert(~isfield(cfg.trialdef, 'offset'));

  if ismember('sample', events.Properties.VariableNames)
    % use the sample number, it is assumed that sample 1 corresponds to the first data sample (and not sample 0)
    begsample = round(events.sample - cfg.trialdef.prestim  * hdr.Fs);
    endsample = round(events.sample + cfg.trialdef.poststim * hdr.Fs);
    offset    = round(-cfg.trialdef.prestim*hdr.Fs) * ones(size(begsample));
  else
    % use the onset, which is in seconds
    begsample = round(events.onset * hdr.Fs + 1 - cfg.trialdef.prestim  * hdr.Fs);
    endsample = round(events.onset * hdr.Fs + 1 + cfg.trialdef.poststim * hdr.Fs);
    offset    = round(-cfg.trialdef.prestim*hdr.Fs) * ones(size(begsample));
  end

elseif isfield(cfg.trialdef, 'duration') && isfield(cfg.trialdef, 'offset')
  % these are mutually exclusive
  assert(~isfield(cfg.trialdef, 'prestim'));
  assert(~isfield(cfg.trialdef, 'poststim'));

  if ismember('sample', events.Properties.VariableNames)
    % use the event sample number, it is assumed that sample 1 corresponds to the first data sample (and not sample 0)
    begsample = round(events.sample);
    endsample = round(events.sample + cfg.trialdef.duration * hdr.Fs);
    offset    = round(cfg.trialdef.offset * hdr.Fs) * ones(size(begsample));
  else
    % use the event onset, which is in seconds
    begsample = round(events.onset * hdr.Fs + 1);
    endsample = round(events.onset * hdr.Fs + cfg.trialdef.duration * hdr.Fs);
    offset    = round(cfg.trialdef.offset * hdr.Fs) * ones(size(begsample));
  end

else
  ft_error('inconsistent specification of the cfg.trialdef options')
end

% construct the minimal required trial definition
trl = table(begsample, endsample, offset);

% remove conflicting columns, as these cannot be used as column headings in the output
fn = intersect({'begsample', 'endsample', 'offset'}, events.Properties.VariableNames);
for i=1:numel(fn)
  ft_warning('ignoring the %s column from the events', fn{i});
  events.(fn{i}) = [];
end

% append the columns with details on the selected events
trl = cat(2, trl, events);

% remove trials that fall outside the data
outside = trl.begsample<1 | trl.endsample>(hdr.nSamples*hdr.nTrials);
if any(outside)
  ft_notice('removing %d trials that extend outside the data', sum(outside));
  trl = trl(~outside, :);
end
