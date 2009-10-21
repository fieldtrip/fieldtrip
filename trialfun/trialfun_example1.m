function [trl, event] = trialfun_example1(cfg)

% TRIALFUN_EXAMPLE1 is an example trial function. It searches for events
% of type "trigger" and specifically for a trigger with value 7, followed
% by a trigger with value 64.

% read the header information and the events from the data
% this should always be done using the generic read_header
% and read_event functions
hdr   = read_fcdc_header(cfg.dataset);
event = read_fcdc_event(cfg.dataset);

% search for "trigger" events
value  = [event(find(strcmp('trigger', {event.type}))).value]';
sample = [event(find(strcmp('trigger', {event.type}))).sample]';

% determine the number of samples before and after the trigger
pretrig  = -cfg.trialdef.pre  * hdr.Fs;
posttrig =  cfg.trialdef.post * hdr.Fs;

% look for the combination of a trigger "7" followed by a trigger "64" 
% for each trigger except the last one
trl = [];
for j = 1:(length(trigger)-1)
  trg1 = trigger(j);
  trg2 = trigger(j+1);
  if trg1==7 && trg2==64
    trlbegin = sample(j) + pretrig;       
    trlend   = sample(j) + posttrig;       
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
  end
end

