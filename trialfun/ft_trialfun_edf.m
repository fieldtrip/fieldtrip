function [trl, event] = ft_trialfun_edf(cfg)

% FT_TRIALFUN_EDF is an example trial function for EDF data. It searches for events
% of type "up" in an analog data channel, as indentified by thresholding. This
% threshold can be a hard threshold, i.e. a numeric, or can flexibly be defined
% depending on the data, for example calculating the 'median' of an analog signal.
%
% You can use this as a template for your own conditial trial definitions.
%
% Use this function by calling 
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset  = string with the filename
%   cfg.trialfun = 'ft_trialfun_edf'
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL

ft_warning('this trial function is only an example, please copy it and adapt it to your specific EDF situation - see http://www.fieldtriptoolbox.org/getting_started/edf');

% read the header information
hdr           = ft_read_header(cfg.dataset);

% read the events from the data
chanindx      = 1; % this should be adapted to your data
detectflank   = 'up';
threshold     = '(3/2)*nanmedian'; % or, e.g., 1/2 times the median for down flanks
event         = ft_read_event(cfg.dataset, 'trigindx', trigindx, 'detectflank', detectflank, 'threshold', threshold);

% define trials around the events
trl           = [];
pretrig       = 1 * hdr.Fs; % e.g., 1 sec before trigger
posttrig      = 2 * hdr.Fs; % e.g., 2 sec after trigger
for i = 1:numel(event)
  offset    = -hdr.nSamplesPre;  % number of samples prior to the trigger
  trlbegin  = event(i).sample - pretrig;
  trlend    = event(i).sample + posttrig;
  newtrl    = [trlbegin trlend offset];
  trl       = [trl; newtrl]; % store in the trl matrix
end
