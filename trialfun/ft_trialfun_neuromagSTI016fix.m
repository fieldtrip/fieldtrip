function [trl, event] = ft_trialfun_neuromagSTI016fix(cfg)

% FT_TRIALFUN_NAME is supposed to fix the error with STI016 in Neuromag
% data. It reads channels with name STI0* (i.e. the channels STI001-STI016)
% and combine the values into a new "STI101" channel. It then use the new
% channel to define trials.
% 
% You would use this function as follows
%   cfg                     = [];   
%   cfg.dataset             = string, containing filename or directory
%   cfg.trialdef.prestim    = pre stim time (in s)
%   cfg.trialdef.poststim   = post stim time (in s)
%   cfg.trialdef.eventvalue = trigger t 
%   cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';
%   cfg                     = definetrial(cfg);
%   data                    = preprocessing(cfg);
%
% You can use this example trial function as template for your own
% conditial trial definitions.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% Get only specific event type
cfg.trialdef.eventtype = ft_getopt(cfg.trialdef, 'eventtype', 'STI101');
cfg.trialdef.eventvalue = ft_getopt(cfg.trialdef, 'eventvalue', []);

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);

% Read trigger channels
chanindx = find(not(cellfun('isempty', strfind(hdr.label,'STI0'))));

event = ft_read_event(cfg.dataset, 'chanindx', chanindx);

% Manually make combined trigger channel
dat = ft_read_data(cfg.dataset , 'chanindx', chanindx); % Read the trigger data
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

