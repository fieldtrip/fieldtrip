function [cfg] = ft_definetrial(cfg)

% FT_DEFINETRIAL defines the segments of data that will be used for
% further processing and analysis, i.e. the pieces of data that will
% be read in by FT_PREPROCESSING. Trials are defined by their begin
% and end sample in the data file and each trial has an offset that
% defines where the relative t=0 point (usually the sample at which
% the trigger is detected) is for that trial.
%
% Use as
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain either
%   cfg.trialdef   = structure with details of trial definition, see below
%   cfg.trialfun   = function name, see below (default = 'trialfun_general')
% and also
%   cfg.dataset    = pathname to dataset from which to read the events
% 
% A call to FT_DEFINETRIAL results in the trial definition "trl" being
% added to the output configuration structure. The trials are defined
% according to the triggers, trials or other events in the data, or
% from a user-specified Matlab function (subsequently referred to as
% the trial function) which returns "trl". The user can specify the
% name of his/her custom trial function that is tailored to the
% experimental paradigm, or use the default trial function
% TRIALFUN_GENERAL.
%
% The trial definition "trl" is an Nx3 matrix, N is the number of trials.
% The first column contains the sample-indices of the begin of each trial 
% relative to the begin of the raw data, the second column contains the 
% sample-indices of the end of each trial, and the third column contains 
% the offset of the trigger with respect to the trial. An offset of 0 
% means that the first sample of the trial corresponds to the trigger. A 
% positive offset indicates that the first sample is later than the trigger, 
% a negative offset indicates that the trial begins before the trigger.
%
% The trial definition "trl" can contain additional columns besides the
% required three that represend begin, end and offset. These additional
% columns can be used by a custom trialfun to provide numeric information
% about each trial such as trigger codes, response latencies, trial
% type and response correctness. The additional columns of the "trl"
% matrix will be represented in data.trialinfo after FT_PREPROCESSING.
%
% Both the default trial function and your own supplied custom trial
% functions in general will call FT_READ_EVENT to collect all event
% information (such as triggers) from your dataset and to select
% pieces of data according to this information.
%
% Simple trial definitions (e.g. based on a single trigger) are supported
% by TRIALFUN_GENERAL, which is the default trial function. This function
% supports the following options
%   cfg.trialdef.eventtype  = 'string'
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
%   cfg.trialdef.prestim    = number, latency in seconds (optional)
%   cfg.trialdef.poststim   = number, latency in seconds (optional)
% If you specify cfg.trialdef.eventtype  = '?' a list with the events in your 
% data file will be displayed on screen.
%
% If you want to read all data from a continous file in a single or in 
% multiple segments, TRIALFUN_GENERAL understands the following options
%    cfg.trialdef.triallength = duration in seconds (can also be 1 or Inf)
%    cfg.trialdef.ntrials     = number of trials (can also be 1 or Inf)
%
% There are also many other complex ways in which you can define data
% pieces of interest, for example based on a conditional sequence of events
% (e.g. stimulus trigger followed by a correct response). For those cases,
% you can supply your own trial function and specify that as
% cfg.trialfun. See belor for pointers to some examples.
% 
% The cfg.trialfun option is a string containing the name of a function
% that you wrote yourself and that FT_DEFINETRIAL will call. The
% function should take the cfg-structure as input and should give a
% NxM matrix with M equal to or larger than 3) in the same format as
% "trl" as the output. You can add extra custom fields to the
% configuration structure to pass as arguments to your own trialfun.
% Furthermore, inside the trialfun you can use the FT_READ_EVENT
% function to get the event information from your data file.
%
% See also FT_PREPROCESSING, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT,
% TRIALFUN_GENERAL, TRIALFUN_EXAMPLE1, TRIALFUN_EXAMPLE2

% Undocumented local options:
% cfg.datafile
% cfg.dataset
% cfg.event
% cfg.trl
% cfg.version

% Copyright (C) 2003, Robert Oostenveld, FCDC
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});

if ~isfield(cfg, 'trl') && (~isfield(cfg, 'trialfun') || isempty(cfg.trialfun))
  % there used to be other system specific trialfuns in previous versions
  % of fieldtrip, but they are deprecated and not included in recent
  % versions any more
  cfg.trialfun = 'trialfun_general';
  warning('no trialfun was specified, using trialfun_general');
end

% create the trial definition for this dataset and condition
if isfield(cfg, 'trl')
  % the trial definition is already part of the configuration
  fprintf('retaining exist trial definition\n');
  trl = cfg.trl;
  if isfield(cfg, 'event')
    fprintf('retaining exist event information\n');
    event = cfg.event;
  else
    event = [];
  end

elseif isfield(cfg, 'trialfun')

  % provide support for xxx and trialfun_xxx when the user specifies cfg.trialfun=xxx
  if isa(cfg.trialfun, 'function_handle') || exist(cfg.trialfun, 'file')
    % evaluate this function, this is the default
  elseif exist(['trialfun_' cfg.trialfun], 'file')
    % prepend trialfun to the function name
    cfg.trialfun = ['trialfun_' cfg.trialfun];
  else
    if ischar(cfg.trialfun)
      error('cannot locate the specified trialfun (%s)', cfg.trialfun)
    else
      error('cannot locate the specified trialfun (%s)', func2str(cfg.trialfun))
    end
  end
  
  % evaluate the user-defined function that gives back the trial definition
  if ischar(cfg.trialfun)
    fprintf('evaluating trialfunction ''%s''\n', cfg.trialfun);
  else
    fprintf('evaluating trialfunction ''%s''\n', func2str(cfg.trialfun));
  end
  
  % determine the number of outpout arguments of the user-supplied trial function
  try
    % the nargout function in Matlab 6.5 and older does not work on function handles
    num = nargout(cfg.trialfun);
  catch
    num = 1;
  end
  if num==1
    % the user-defined function only gives back the trial definition
    trl   = feval(cfg.trialfun, cfg);
    event = [];
  else
    % the user-defined function also gives back detailed information about
    % conditions, reaction time or any other information
    [trl, event] = feval(cfg.trialfun, cfg);
  end
else
  error('no trialfunction specified, see FT_DEFINETRIAL for help');
end

if isfield(cfg, 'trialdef') && isfield(cfg.trialdef, 'eventtype') && isequal(cfg.trialdef.eventtype, '?')
  % give a gentle message instead of an error
  fprintf('no trials have been defined yet, see FT_DEFINETRIAL for further help\n');
elseif size(trl,1)<1
  error('no trials were defined, see FT_DEFINETRIAL for help');
end

% add the new trials and events to the output configuration
fprintf('found %d events\n', length(event));
cfg.event = event;
fprintf('created %d trials\n', size(trl,1));
cfg.trl = trl;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble callinfo

