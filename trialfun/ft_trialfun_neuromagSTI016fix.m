function [trl, event] = ft_trialfun_neuromagSTI016fix(cfg)

% FT_TRIALFUN_NEUROMAGSTI106FIX is supposed to fix the error with STI016 in
% Neuromag/Elekta/MEGIN data. It reads the channels STI001 up to STI016, combines the
% values into a new "STI101" channel and then uses the new channel to define trials.
%
% Use this function by calling
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset             = string, containing filename or directory
%   cfg.trialdef.prestim    = pre stimulus time in s
%   cfg.trialdef.poststim   = post stimulus time in seconds
%   cfg.trialdef.eventvalue = list with trigger values
%   cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL

% Undocumented option:
%  cfg.checkmaxfilter      = check that MaxFilter has run (default = yes).

% Get only specific event type
cfg.trialdef.eventtype = ft_getopt(cfg.trialdef, 'eventtype', 'STI101');
cfg.trialdef.eventvalue = ft_getopt(cfg.trialdef, 'eventvalue', []);

% MaxFilter option
cfg.checkmaxfilter = ft_getopt(cfg.trialdef, 'checkmaxfilter');

% read the header information, needed for the sampling rate and channel labels
hdr = ft_read_header(cfg.dataset, 'checkmaxfilter', cfg.checkmaxfilter);

% Read trigger channels
chanindx = find(not(cellfun('isempty', strfind(hdr.label,'STI0'))));
event = ft_read_event(cfg.dataset, 'chanindx', chanindx, 'checkmaxfilter', cfg.checkmaxfilter);
    
% Manually make combined trigger channel
dat = ft_read_data(cfg.dataset , 'chanindx', chanindx, 'checkmaxfilter', cfg.checkmaxfilter); % Read the trigger data
allsti = dat(1:16,:)==5;
trig = allsti(1,:)*2^0 + ...
  allsti(2,:)*2^1 + ...
  allsti(3,:)*2^2 + ...
  allsti(4,:)*2^3 + ...
  allsti(5,:)*2^4 + ...
  allsti(6,:)*2^5 + ...
  allsti(7,:)*2^6 + ...
  allsti(8,:)*2^7 + ...
  allsti(9,:)*2^8 + ...
  allsti(10,:)*2^9 + ...
  allsti(11,:)*2^10 + ...
  allsti(12,:)*2^11 + ...
  allsti(13,:)*2^12 + ...
  allsti(14,:)*2^13 + ...
  allsti(15,:)*2^14 + ...
  allsti(16,:)*2^15;

% event = struct();
for j=find(diff([0 trig])>0)
  event(end+1).type   = 'STI101';
  event(end  ).sample = j;          % assign the sample at which the trigger has gone up
  event(end  ).value  = trig(j);    % assign the trigger value just _after_ going up
end

%
value  = [event.value]';
sample = [event.sample]';

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig =  round(cfg.trialdef.poststim * hdr.Fs);

if isempty(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = unique(value);
end

% look for the triggers
trl = [];
for j = 1:length(value)
  if strcmp(cfg.trialdef.eventtype, event(j).type)
    trg = value(j);
    if any(cfg.trialdef.eventvalue==trg)
      trlbegin = sample(j) + pretrig;
      trlend   = sample(j) + posttrig;
      offset   = pretrig;
      newtrl   = [trlbegin trlend offset, trg];
      trl      = [trl; newtrl];
    end
  end
end
