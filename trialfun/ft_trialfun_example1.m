function [trl, event] = ft_trialfun_example1(cfg)

% FT_TRIALFUN_EXAMPLE1 is an example trial function. It searches for events
% of type "trigger" and specifically for a trigger with value 7, followed
% by a trigger with value 64.
% 
% You can use this example trial function as template for your own
% conditial trial definitions.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% read the header and determine the channel number corresponding with the EMG
hdr         = ft_read_header(cfg.headerfile);
chanindx    = strmatch('EMGlft', hdr.label);
 
if length(chanindx)>1
  error('only one EMG channel supported');
end
 
% read all data of the EMG channel, assume continuous file format
begsample = 1;
endsample = hdr.nSamples*hdr.nTrials;
emg       = ft_read_data(cfg.datafile, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);
 
% apply filtering, hilbert transformation and boxcar convolution (for smoothing)
emgflt      = ft_preproc_highpassfilter(emg, hdr.Fs, 10);     % highpassfilter (helper function from fieldtrip)
emghlb      = abs(hilbert(emgflt')');                     % hilbert transform
emgcnv      = conv2([1], ones(1,hdr.Fs), emghlb, 'same'); % smooth using convolution
emgstd      = ft_preproc_standardize(emgcnv);                % z-transform (helper function from fieldtrip)
emgtrl      = emgstd>0;                                   % detect the muscle activity by thresholding
emgtrl      = diff(emgtrl, [], 2);
 
% find the onset and offset
emgon       = find(emgtrl(:)== 1);
emgoff      = find(emgtrl(:)==-1);
 
trl(:,1) = emgon (:) + hdr.Fs*0.5;  % as a consequence of the convolution with a one-second boxcar
trl(:,2) = emgoff(:) - hdr.Fs*0.5;  % as a consequence of the convolution with a one-second boxcar
trl(:,3) = 0;

% read the header information and the events from the data
% this should always be done using the generic read_header
% and read_event functions
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% search for "trigger" events
trigger = [event(strcmp('trigger', {event.type})).value]';
sample  = [event(strcmp('trigger', {event.type})).sample]';

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
posttrig =  round(cfg.trialdef.post * hdr.Fs);

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

