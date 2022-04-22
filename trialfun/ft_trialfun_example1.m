function [trl, event] = ft_trialfun_example1(cfg)

% FT_TRIALFUN_EXAMPLE1 is an example trial function. It searches for events
% of type "trigger" and specifically for a trigger with value 7, followed
% by a trigger with value 64.
%
% You can use this as a template for your own conditial trial definitions.
%
% Use this function by calling
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset           = string with the filename
%   cfg.trialfun          = 'ft_trialfun_example1'
%   cfg.trialdef.prestim  = number, in seconds
%   cfg.trialdef.poststim = number, in seconds
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% search for "trigger" events
value  = [event(find(strcmp('trigger', {event.type}))).value]';
sample = [event(find(strcmp('trigger', {event.type}))).sample]';

% determine the number of samples before and after the trigger
prestim  = -round(cfg.trialdef.prestim  * hdr.Fs);
poststim =  round(cfg.trialdef.poststim * hdr.Fs);

% look for the combination of a trigger "7" followed by a trigger "64"
% for each trigger except the last one
trl = [];
for j = 1:(length(value)-1)
  trg1 = value(j);
  trg2 = value(j+1);
  if trg1==7 && trg2==64
    trlbegin = sample(j) + prestim;
    trlend   = sample(j) + poststim;
    offset   = prestim;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
  end
end
