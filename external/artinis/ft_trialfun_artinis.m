function [trl, event] = ft_trialfun_artinis(cfg)

% FT_TRIALFUN_ARTINIS is adjoining the result of ft_trialfun_general and 
% those events found by FT_TRIALFUN_GENERAL.
%
% The trialdef structure can contain the following specifications
%   cfg.trialdef.eventtype  = 'string'
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
%   cfg.trialdef.oxyproj    = 'string', indicating an oxyproj-file, in
%                             which information about the events for this
%                             oxy3-file are stored
%   cfg.trialdef.prestim    = latency in seconds (optional)
%   cfg.trialdef.poststim   = latency in seconds (optional)
%
% If you want to read all data from a continous file in segments, you can specify
%    cfg.trialdef.triallength = duration in seconds (can be Inf)
%    cfg.trialdef.ntrials     = number of trials
%
% If you specify
%   cfg.trialdef.eventtype  = '?'
% a list with the events in your datafile will be displayed on screen.
%
% If you specify
%   cfg.trialdef.eventtype = 'gui'
% a graphical user interface will allow you to select events of interest.
%
% See also FT_TRIALFUN_GENERAL, FT_DEFINETRIAL, FT_PREPROCESSING

% You are using the FieldTrip NIRS toolbox developed and maintained by 
% Artinis Medical Systems (http://www.artinis.com). For more information
% on FieldTrip, see http://www.fieldtriptoolbox.org
% 
% This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 
% International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% 
% Creative Commons Attribution-ShareAlike 4.0 International License:
% -----------------------------------
% You are free to:
% 
%     Share — copy and redistribute the material in any medium or format
%     Adapt — remix, transform, and build upon the material
%     for any purpose, even commercially.
% 
%     The licensor cannot revoke these freedoms as long as you follow the 
%     license terms.
% 
% Under the following terms:
% 
%     Attribution — You must give appropriate credit, provide a link to 
%                    the license, and indicate if changes were made. You 
%                    may do so in any reasonable manner, but not in any way 
%                    that suggests the licensor endorses you or your use.
% 
%     ShareAlike — If you remix, transform, or build upon the material, 
%                   you must distribute your contributions under the same 
%                   license as the original.
% 
%     No additional restrictions — You may not apply legal terms or 
%                                   technological measures that legally 
%                                   restrict others from doing anything the 
%                                   license permits.
% 
% -----------------------------------
% 
% This toolbox is not to be used for medical or clinical purposes.
% 
% Copyright (c) 2016 by Artinis Medical Systems.
% Contact: askforinfo@artinis.com
%
% Main programmer: 
% Jörn M. Horschig, Artinis Medical Systems BV, http://www.artinis.com
% $Id$

if ~isfield(cfg.trialdef, 'oxyproj') || ~strcmp(cfg.trialdef.oxyproj(end-7:end), '.oxyproj') || ~exist(cfg.trialdef.oxyproj, 'file') 
  error('You need to provide an existing oxyproj-file as input to use this function (cfg.trialdef.oxyproj)')
end

cfg.eventformat = ft_getopt(cfg, 'eventformat');
cfg.headerformat = ft_getopt(cfg, 'headerformat');
cfg.dataformat = ft_getopt(cfg, 'dataformat');

fprintf('reading the events from ''%s''\n', cfg.headerfile);
event = ft_read_event(cfg.headerfile, 'headerformat', cfg.headerformat, 'eventformat', cfg.eventformat, 'dataformat', cfg.dataformat);

fprintf('reading the events from ''%s''\n', cfg.trialdef.oxyproj);
measurement = parseoxyproj(cfg.trialdef.oxyproj, cfg.datafile);

if numel(measurement)>1
  warning('More than one measurement found in the oxyproj file - removing all but the first one\n');
  measurement(2:end) = [];
end

% read the header for obtaining the sample rate
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);

% create the trl and event-structure with the new information
trl = [];
for i=1:numel(measurement.events.names)
  add2trl = strcmp(cfg.trialdef.eventtype, 'oxyproj') && strcmp(cfg.trialdef.eventvalue, measurement.events.names{i});
  for e=1:numel(measurement.events.onsets{i})    
    % add event
    event(end+1).type      = 'oxyproj';
    event(end).sample    = measurement.events.onsets{i}(e);
    event(end).value     = measurement.events.names{i};
    event(end).offset    = 0;
    event(end).duration  = measurement.events.durations{i}(e);    
    if add2trl
      trloff = round(-cfg.trialdef.prestim * hdr.Fs);
      trlbeg = event(end).sample+trloff;
      % this will not work if prestim was not defined, the code will then crash
      trldur = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
      trlend = trlbeg + trldur;      
      trl(end+1, :) = [trlbeg trlend trloff];
    end
  end
end

% create trials
cfg.event = event;
[trl, event] = ft_trialfun_general(cfg);
