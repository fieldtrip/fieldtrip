function [trl] = ft_trialfun_emgdetect(cfg)

% Note that there are some parameters, like the EMG channel name and the
% processing that is done on the EMG channel data, which are hardcoded in
% this trial function. You should change these parameters if necessary.

% read the header and determine the channel number corresponding with the EMG
hdr         = ft_read_header(cfg.headerfile);
chanindx    = strmatch('EMGlft', hdr.label);

if length(chanindx)>1
  error('only one EMG channel supported');
end

% read all data of the EMG channel, assume continuous file format
emg = ft_read_data(cfg.datafile, 'header', hdr, ...
              'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials, ...
              'chanindx', chanindx, 'checkboundary', false);

% apply filtering, hilbert transformation and boxcar convolution (for smoothing)
emgflt      = ft_preproc_highpassfilter(emg, hdr.Fs, 10); % highpassfilter
emghlb      = abs(hilbert(emgflt')');                     % hilbert transform
emgcnv      = conv2([1], ones(1,hdr.Fs), emghlb, 'same'); % smooth using convolution
emgstd      = ft_preproc_standardize(emgcnv, 2);          % z-transform, i.e. mean=0 and stdev=1
emgtrl      = emgstd>0;                                   % detect the muscle activity
emgtrl      = diff(emgtrl, [], 2);

emgon       = find(emgtrl(:)== 1);
emgoff      = find(emgtrl(:)==-1);

trl(:,1) = emgon (:) + hdr.Fs*0.5;  % as a consequence of the convolution with a one-second boxcar
trl(:,2) = emgoff(:) - hdr.Fs*0.5;  % as a consequence of the convolution with a one-second boxcar
trl(:,3) = 0;
