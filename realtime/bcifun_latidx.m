function cmd = bcifun_latidx(cfg,data)

% BCIFUN_LATIDX extracts power, computes index and extracts commands
%
% this function outputs the lateralization index:
%
% log A / B
%
% for channel averages A and B
%
% Options:
% cfg.foilim  = the frequency band of interest (default = [8 14])
% cfg.avgchan = a cell array that specifies over which channels to average;
%               this can be either [], meaning no averaging or a cell-array
%               of channel index vectors; default behaviour is to average
%               over the first n/2 channels and the second n/2 channels;
%
% Copyright (C) 2009, Marcel van Gerven
%
% $Log: bcifun_latidx.m,v $
% Revision 1.2  2009/06/30 12:38:33  roboos
% *** empty log message ***
%
% Revision 1.1  2009/04/23 12:50:49  marvger
% update of the BCI realtime code
%

% set the default configuration options
if ~isfield(cfg, 'foilim'),  cfg.foilim = [8 14];     end
if ~isfield(cfg, 'avgchan')
  n = length(data.label);
  cfg.avgchan = { 1:floor(n/2) (floor(n/2)+1):n };
end

cmd = nan(1,length(data.trial));

for trllop=1:length(data.trial)

  % set up the spectral estimator
  specest = spectrum.welch('Hamming', min(data.fsample, size(data.trial{trllop},2)));

  nchan = length(data.label);

  if nchan < 2 || length(cfg.avgchan) ~= 2
    error('neurofeedback application expects exactly two channel groups');
  end

  features = zeros(1,nchan);
  for i=1:nchan
    est = psd(specest, data.trial{trllop}(i,:), 'Fs', data.fsample);
    features(i) = mean(est.Data(est.Frequencies >= cfg.foilim(1) & est.Frequencies <= cfg.foilim(2)));
  end

  % average over channel subsets
  if ~isempty(cfg.avgchan)
    cfeatures  = zeros(1,length(cfg.avgchan));
    for j=1:length(cfg.avgchan)
      cfeatures(j) = mean(features(cfg.avgchan{j}));
    end
    features = cfeatures;
  end

  cmd(trllop) = log10(features(1) / features(2));

end % looping over trials
