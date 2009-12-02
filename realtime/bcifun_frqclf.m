function cmd = bcifun_frqclf(cfg,data)

% BCIFUN_FREQCLF classifies data based on powerspectra. It is assumed that cfg.count represents
% the index of the current example. The bcifun makes use of persistent
% states.
%
% Options:
% cfg.foilim  = the frequency band of interest (default = [8 14])
% cfg.clfproc = the classification procedure (default = {standardizer() svmmethod()})
% cfg.ntrain  = the number of training examples (default = 40)
%
% Copyright (C) 2009, Marcel van Gerven
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

persistent DATA;
persistent LABELS;
persistent CLFPROC;

% set the default configuration options
if ~isfield(cfg, 'foilim'),  cfg.foilim = [8 14];     end
if ~isfield(cfg, 'clfproc'), cfg.clfproc = clfproc({standardizer() svmmethod()}); end
if ~isfield(cfg, 'ntrain'),  cfg.ntrain = 20; end

cmd = nan(1,length(data.trial));

for trllop=1:length(data.trial)

  nchan = length(data.label);

  % set up the spectral estimator; FIXME: this code should be replaced by
  % native FieldTrip code
  specest = spectrum.welch('Hamming', min(data.fsample, size(data.trial{trllop},2)));
  features = zeros(1,nchan);
  for i=1:nchan
    est = psd(specest, data.trial{trllop}(i,:), 'Fs', data.fsample);
    features(i) = mean(est.Data(est.Frequencies >= cfg.foilim(1) & est.Frequencies <= cfg.foilim(2)));
  end

  % initialization
  if cfg.count == 1
    DATA    = zeros(cfg.ntrain,length(features));
    LABELS  = zeros(cfg.ntrain,1);
    CLFPROC = cfg.clfproc;
  end

  if cfg.count <= cfg.ntrain % training mode

    % save features and class label
    DATA(cfg.count,:) = features;
    LABELS(cfg.count) = data.cfg.trl(1,4);

    if cfg.count == cfg.ntrain % train classifier
      fprintf('training classification procedure\n');
      CLFPROC = CLFPROC.train(DATA,LABELS);
    end

    cmd(trllop) = nan;

  else % test mode

    cmd(trllop) = CLFPROC.predict(features);

  end

end % looping over trials
