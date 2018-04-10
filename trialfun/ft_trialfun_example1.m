function [trl, event] = ft_trialfun_example1(cfg)

% FT_TRIALFUN_EXAMPLE1 is an example trial function. It searches for events
% of type "trigger" and specifically for a trigger with value 7, followed
% by a trigger with value 64.
% 
% You would use this function as follows
%   cfg           = [];   
%   cfg.dataset   = string, containing filename or directory
%   cfg.trialfun  = 'ft_trialfun_example1';
%   cfg           = definetrial(cfg);
%   data          = preprocessing(cfg);
%
% You can use this example trial function as template for your own
% conditial trial definitions.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% search for "trigger" events
value  = [event(find(strcmp('trigger', {event.type}))).value]';
sample = [event(find(strcmp('trigger', {event.type}))).sample]';

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
posttrig =  round(cfg.trialdef.post * hdr.Fs);

% look for the combination of a trigger "7" followed by a trigger "64" 
% for each trigger except the last one
trl = [];
for j = 1:(length(value)-1)
  trg1 = value(j);
  trg2 = value(j+1);
  if trg1==7 && trg2==64
    trlbegin = sample(j) + pretrig;       
    trlend   = sample(j) + posttrig;       
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
  end
end

