function [trl, event] = ft_trialfun_imotions(cfg)

% FT_TRIALFUN_IMOTIONS makes a trial definition for an iMotions event structure
%
% The trialdef structure can contain the following specifications
% cfg.trialdef.eventtype  = string or cell-array of strings (default = 'StimulusName')
% cfg.trialdef.eventvalue = string or cell-array of strings (default = [])
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL

assert(isfield(cfg, 'event'), 'the configuration should contain the event structure');
assert(isfield(cfg, 'fsample'), 'the configuration should contain the sampling frequency');

% set the defaults
cfg.trialdef            = ft_getopt(cfg, 'trialdef', []);
cfg.trialdef.eventtype  = ft_getopt(cfg.trialdef, 'eventtype', 'StimulusName');
cfg.trialdef.eventvalue = ft_getopt(cfg.trialdef, 'eventvalue', []);
cfg.trialdef.offset     = ft_getopt(cfg.trialdef, 'offset', 'absolute'); % absolute or relative

event = cfg.event;

% start by selecting all events
sel = true(size(event));
% make a subselection of events
if ~isempty(cfg.trialdef.eventtype)
  sel = sel & strcmp({event.type}, cfg.trialdef.eventtype);
end
if ~isempty(cfg.trialdef.eventvalue)
  sel = sel & strcmp({event.value}, cfg.trialdef.eventvalue);
end

sample   = [event(sel).sample];
duration = [event(sel).duration];
type     = {event(sel).type};
value    = {event(sel).value};

begsample = sample(:);
endsample = sample(:)+duration(:)-1;

switch cfg.trialdef.offset
  case 'relative'
    % start of each trial/segment is t=0
    offset    = zeros(size(begsample));
  case 'absolute'
    % start of recording is t=0
    offset    = begsample - 1;
end

type  = type(:);
value = value(:);
trl   = table(begsample, endsample, offset, type, value);



