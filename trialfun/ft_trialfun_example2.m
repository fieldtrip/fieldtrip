function [dat] = ft_trialfun_example2(cfg)
 
% FT_TRIALFUN_EXAMPLE2 is an example trial function that detects muscle
% activity in an EMG channel and defines variable length trials from the
% EMG onset up to the EMG offset.
%
% You would use this function as follows
%   cfg           = [];   
%   cfg.dataset   = string, containing filename or directory
%   cfg.trialfun  = 'ft_trialfun_example2';
%   cfg           = definetrial(cfg);
%   data          = preprocessing(cfg);
%
% Note that there are some parameters, like the EMG channel name and the
% processing that is done on the EMG channel data, which are hardcoded in
% this trial function. You should change these parameters if necessary.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% read the header and determine the channel number corresponding with the EMG
hdr         = ft_read_header(cfg.headerfile);
chanindx    = strmatch('EMGlft', hdr.label);
 
if length(chanindx)>1
  ft_error('only one EMG channel supported');
end
 
% read all data of the EMG channel, assume continuous file format
begsample = 1;
endsample = hdr.nSamples*hdr.nTrials;
emg       = ft_read_data(cfg.datafile, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);
 
% apply filtering, hilbert transformation and boxcar convolution (for smoothing)
emgflt      = ft_preproc_highpassfilter(emg, hdr.Fs, 10); % highpassfilter (helper function from fieldtrip/preproc)
emghlb      = abs(hilbert(emgflt')');                     % hilbert transform
emgcnv      = conv2([1], ones(1,hdr.Fs), emghlb, 'same'); % smooth using convolution
emgstd      = ft_preproc_standardize(emgcnv);             % z-transform (helper function from fieldtrip/preproc)
emgtrl      = emgstd>0;                                   % detect the muscle activity by thresholding
emgtrl      = diff(emgtrl, [], 2);
 
% find the onset and offset
emgon       = find(emgtrl(:)== 1);
emgoff      = find(emgtrl(:)==-1);
 
trl(:,1) = emgon (:) + hdr.Fs*0.5;  % as a consequence of the convolution with a one-second boxcar
trl(:,2) = emgoff(:) - hdr.Fs*0.5;  % as a consequence of the convolution with a one-second boxcar
trl(:,3) = 0;

