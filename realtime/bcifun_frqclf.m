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
% $Log: bcifun_frqclf.m,v $
% Revision 1.6  2009/07/01 09:48:52  marvger
% changed allcmd back to cmd
%
% Revision 1.5  2009/06/30 12:38:50  roboos
% made consistent with other fieldtrip behaviour
%
% Revision 1.4  2009/06/14 12:24:08  marvger
% changed classify to predict
%
% Revision 1.3  2009/04/29 10:27:54  marvger
% update
%
% Revision 1.2  2009/04/23 12:53:27  marvger
% changed CLFPROC.test to CLFPROC.classify
%
% Revision 1.1  2009/04/23 12:50:49  marvger
% update of the BCI realtime code
%

persistent DATA;
persistent LABELS;
persistent CLFPROC;

% set the default configuration options
if ~isfield(cfg, 'foilim'),  cfg.foilim = [8 14];     end
if ~isfield(cfg, 'clfproc'), cfg.clfproc = clfproc({standardizer() svmmethod()}); end
if ~isfield(cfg, 'ntrain'),  cfg.ntrain = 20; end

cmd = nan(1,length(data.trial));

for trllop=1:length(data.trial)
  
  % set up the spectral estimator
  specest = spectrum.welch('Hamming', min(data.fsample, size(data.trial{trllop},2)));

  nchan = length(data.label);
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
