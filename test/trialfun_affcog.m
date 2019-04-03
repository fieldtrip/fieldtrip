
function [trl, event] = trialfun_affcog(cfg)

%% the first part is common to all trial functions
% read the header (needed for the samping rate) and the events
hdr        = ft_read_header(cfg.headerfile);
event      = ft_read_event(cfg.headerfile);

%% from here on it becomes specific to the experiment and the data format
% for the events of interest, find the sample numbers (these are integers)
% for the events of interest, find the trigger values (these are strings in the case of BrainVision)
EVsample   = [event.sample]';
EVvalue    = {event.value}';

% select the target stimuli
Word = find(strcmp('S141', EVvalue)==1);

% for each word find the condition
for w = 1:length(Word)
  % code for the judgement task: 1 => Affective; 2 => Ontological;
  if strcmp('S131', EVvalue{Word(w)+1}) == 1
    task(w,1) = 1;
  elseif strcmp('S132', EVvalue{Word(w)+1}) == 1
    task(w,1) = 2;
  end
end

PreTrig   = round(0.2 * hdr.Fs);
PostTrig  = round(1 * hdr.Fs);

begsample = EVsample(Word) - PreTrig;
endsample = EVsample(Word) + PostTrig;

offset = -PreTrig*ones(size(endsample));

%% the last part is again common to all trial functions
% return the trl matrix (required) and the event structure (optional)
trl = [begsample endsample offset task];

end % function
