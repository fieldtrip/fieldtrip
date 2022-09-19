function [stat, cfg] = ft_statistics_mvpa(cfg, dat, design)

% FT_STATISTICS_MVPA performs multivariate pattern classification or regression using
% the MVPA-Light toolbox. The function supports cross-validation, searchlight
% analysis, generalization, nested preprocessing, a variety of classification and
% regression metrics, as well as statistical testing of these metrics. This function
% should not be called directly, instead you should call the function that is
% associated with the type of data on which you want to perform the test.
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
%
% where the data is obtained from FT_TIMELOCKANALYSIS, FT_FREQANALYSIS or
% FT_SOURCEANALYSIS respectively, or from FT_TIMELOCKGRANDAVERAGE,
% FT_FREQGRANDAVERAGE or FT_SOURCEGRANDAVERAGE respectively
% and with cfg.method = 'mvpa'
%
% The configuration options that can be specified are:
%   cfg.features        = specifies the name or index of the dimension(s)
%                         that serve(s) as features for the classifier or
%                         regression model. Dimensions that are not
%                         samples or features act as search
%                         dimensions. For instance, assume the data is a
%                         3D array of size [samples x channels x time].
%                         If mvpa.features = 2, the channels serve as
%                         features. A classification is then performed for
%                         each time point (we call time a searchlight
%                         dimension). Conversely, if mvpa.features = 3, the
%                         time points serve as features. A classification
%                         is performed for each channel (channel is a
%                         searchlight dimension).
%                         If mvpa.features = [], then all non-sample
%                         dimensions serve as searchlight dimensions.
%                         If the dimensions have names (ie cfg.dimord
%                         exists), then instead of numbers the feature can
%                         be specified as a string (e.g. 'chan').
%                         (default 2)
%   cfg.generalize      = specifies the name or index of the dimensions
%                         that serves for generalization (if any). For
%                         instance, if the data is [samples x channels x
%                         time], and mvpa.generalize = 3, a time x time
%                         generalization is performed. If mvpa.generalize =
%                         2, a electrode x electrode generalization is
%                         performed. mvpa.generalize must refer to a
%                         searchlight dimension, therefore its value must
%                         be different from the value of mvpa.features.
%                         (default [])
%
% The configuration contains a substruct cfg.mvpa that contains detailed
% options for the MVPA. Possible fields
%   cfg.mvpa.classifier  = string specifying the classifier
%                         Available classifiers:
%                         'ensemble'     Ensemble of classifiers. Any of the other
%                                        classifiers can be used as a learner.
%                         'kernel_fda'   Kernel Fisher Discriminant Analysis
%                         'lda'          Regularized linear discriminant analysis
%                                        (LDA) (for two classes)
%                         'logreg'       Logistic regression
%                         'multiclass_lda' LDA for more than two classes
%                         'naive_bayes'  Naive Bayes
%                         'svm'          Support Vector Machine (SVM)
%                         More details on the classifiers: https://github.com/treder/MVPA-Light#classifiers-for-two-classes-
%                         Additionally, you can choose 'libsvm' or
%                         'liblinear' as a model. They provide interfaces
%                         for logistic regression, SVM, and Support Vector
%                         Regression. Note that they can act as either
%                         classifiers or regression models. An installation
%                         of LIBSVM or LIBLINEAR is required.
%   cfg.mvpa.model       = string specifying the regression model. If a
%                         regression model has been specified,
%                         cfg.mvpa.classifier should be empty (and vice
%                         versa). If neither a classifier nor regression
%                         model is specified, a LDA classifier is used by
%                         default.
%
%                         Available regression models:
%                         'ridge         Ridge regression
%                         'kernel_ridge' Kernel Ridge regression
%                         More details on the regression models: https://github.com/treder/MVPA-Light#regression-models-
%   cfg.mvpa.metric      = string, classification or regression metric, or
%                         cell array with multiple metrics.
%                         Classification metrics: accuracy auc confusion
%                             dval f1 kappa precision recall tval
%                         Regression metrics: mae mse r_squared
%
%   cfg.mvpa.hyperparameter = struct, structure with hyperparameters for the
%                         classifier or regression model (see HYPERPARAMETERS below)
%   cfg.mvpa.feedback       = 'yes' or 'no', whether or not to print feedback on the console (default 'yes')
%
% To obtain a realistic estimate of classification performance, cross-validation
% is used. It is controlled by the following parameters:
%   cfg.mvpa.cv          = string, cross-validation type, either 'kfold', 'leaveout'
%                         'holdout', or 'predefined'. If 'none', no cross-validation is
%                         used and the model is tested on the training
%                         set. (default 'kfold')
%   cfg.mvpa.k           = number of folds in k-fold cross-validation (default 5)
%   cfg.mvpa.repeat      = number of times the cross-validation is repeated
%                         with new randomly assigned folds (default 5)
%   cfg.mvpa.p           = if cfg.cv is 'holdout', p is the fraction of test
%                         samples (default 0.1)
%   cfg.mvpa.stratify    = if 1, the class proportions are approximately
%                         preserved in each test fold (default 1)
%   cfg.mvpa.fold        = if cv='predefined', fold is a vector of length
%                         #samples that specifies the fold each sample belongs to
%
% More information about each classifier is found in the documentation of
% MVPA-Light (github.com/treder/MVPA-Light/).
%
% HYPERPARAMETERS:
% Each classifier comes with its own set of hyperparameters, such as
% regularization parameters and the kernel. Hyperparameters can be set
% using the cfg.mvpa.hyperparameter substruct. For instance, in LDA,
% cfg.mvpa.hyperparameter = 'auto' sets the lambda regularization parameter.
%
% The specification of the hyperparameters is found in the training function
% for each model at github.com/treder/MVPA-Light/tree/master/model
% If a hyperparameter is not specified, default values are used.
%
% SEARCHLIGHT ANALYSIS:
% Data dimensions that are not samples or features serve as 'search
% dimensions'. For instance, if the data is [samples x chan x time]
% and mvpa.features = 'time', then the channel dimension serves as search
% dimension: a separate analysis is carried out for each channel. Instead
% of considering each channel individually, a searchlight can be defined
% such that each channel is used together with its neighbours. Neighbours
% can be specified using the cfg.neighbours field:
%
%   cfg.neighbours   = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%   cfg.timwin       = integer, if MVPA is performed for each time point,
%                      timwin specfies the total size of the time window
%                      that is considered as features.
%                      Example: for cfg.timwin = 3 a given time point is considered
%                      together with the immediately preceding and following
%                      time points. Increasing timwin typially
%                      leads to smoother results along the time axis.
%   cfg.freqwin      = integer, acts like cfg.timwin but across frequencies
%
% This returns:
%   stat.metric = this contains the requested metric
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS, FT_SOURCESTATISTICS,
% FT_STATISTICS_ANALYTIC, FT_STATISTICS_STATS, FT_STATISTICS_MONTECARLO, FT_STATISTICS_CROSSVALIDATE

% Copyright (C) 2019-2021, Matthias Treder and Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

ft_hastoolbox('mvpa-light', 1);

% do a sanity check on the input data
assert(isnumeric(dat),    'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');
assert(isnumeric(design), 'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');

%% cfg: set defaults
cfg.generalize        = ft_getopt(cfg, 'generalize',   []);
cfg.timwin            = ft_getopt(cfg, 'timwin',       []);
cfg.freqwin           = ft_getopt(cfg, 'freqwin',      []);
cfg.neighbours        = ft_getopt(cfg, 'neighbours',   []);
cfg.connectivity      = ft_getopt(cfg, 'connectivity', []); % the default is dealt with below

cfg.mvpa              = ft_getopt(cfg, 'mvpa',       []);
cfg.mvpa.model        = ft_getopt(cfg.mvpa, 'model', []);
if isempty(cfg.mvpa.model)
  cfg.mvpa.classifier = ft_getopt(cfg.mvpa, 'classifier', 'lda');
  if strcmp(cfg.mvpa.classifier, 'naive_bayes')
    cfg.mvpa.append   = ft_getopt(cfg.mvpa, 'append', 1);
  end
  cfg.mvpa.metric     = ft_getopt(cfg.mvpa, 'metric', 'accuracy');
else
  cfg.mvpa.metric     = ft_getopt(cfg.mvpa, 'metric', 'mae');
end
cfg.mvpa.neighbours   = ft_getopt(cfg.mvpa, 'neighbours',  []);
cfg.mvpa.feedback     = ft_getopt(cfg.mvpa, 'feedback',   'yes');

has_dimord = isfield(cfg, 'dimord');
has_dim    = isfield(cfg, 'dim');
if ~has_dimord || ~has_dim
  ft_warning('fields dim or dimord are missing, MVPA may not work correctly')
end

% flip dimensions such that the number of trials comes first
dat = dat.';

if has_dim
  % reshape because MVPA-Light expects the original multi-dimensional array
  dat = reshape(dat, [size(dat,1) cfg.dim]);
end

%% defaults for cfg.features
if ~isfield(cfg, 'features')
  % no features provided, we have to guess
  if ~isempty(cfg.neighbours) || ~isempty(cfg.connectivity)
    if isempty(cfg.timwin)
      cfg.features = 3; % dimension 3 = time (usually)
    end
  elseif ~isempty(cfg.timwin) || ~isempty(cfg.freqwin)
    cfg.features = 2; % dimension 2 = chan (usually)
  end
end
cfg.features = ft_getopt(cfg, 'features', []);

%% backward compatibility
cfg.mvpa = ft_checkconfig(cfg.mvpa, 'renamed', {'param', 'hyperparameter'});
ft_checkconfig(cfg, 'deprecated', {'timextime' 'searchlight'});
ft_checkconfig(cfg.mvpa, 'deprecated', {'balance' 'normalise' 'replace'});

if isfield(cfg,'timextime') && strcmp(cfg.timextime, 'yes')
  cfg.generalize = 'time';
end
if isfield(cfg,'searchlight') && strcmp(cfg.searchlight, 'yes')
  cfg.features = 3;
end
if isfield(cfg,'normalise') && ~isfield(cfg.mvpa,'preprocess')
  cfg = add_to_preprocess(cfg, cfg.normalise);
end
if isfield(cfg,'balance') && ~isempty(cfg.balance)
  cfg = add_to_preprocess(cfg, cfg.balance);
end

%% get dimension names
if has_dimord
  dimtok = tokenize(cfg.dimord, '_');
  cfg.mvpa.dimension_names = ft_getopt(cfg.mvpa, 'dimension_names', [{'samples'} dimtok]);
end

%% convert features and generalize from char to integers
if ischar(cfg.features) || iscell(cfg.features)
  if ~iscell(cfg.features),  feat = {cfg.features};
  else, feat = cfg.features;
  end
  assert(has_dimord, 'if cfg.features is a string then cfg.dimord must exist')
  cfg.features = zeros(1, numel(feat));
  for ix = 1:numel(feat)
    find_ix = find(ismember(cfg.mvpa.dimension_names, feat{ix}));
    assert(~isempty(find_ix), sprintf('''%s'' specified as feature but it is not found in cfg.dimord', feat{ix}))
    cfg.features(ix) = find_ix;
  end
end

if ischar(cfg.generalize)
  assert(has_dimord, 'if cfg.generalize is a string then cfg.dimord must exist')
  cfg.generalize = find(ismember(cfg.mvpa.dimension_names, cfg.generalize));
  if isempty(cfg.generalize)
    ft_error(sprintf('cfg.generalize = ''%s'' is not contained in cfg.dimord', cfg.generalize))
  end
end

cfg.mvpa.feature_dimension          = cfg.features;
cfg.mvpa.generalization_dimension   = cfg.generalize;

% names of search dimensions
dimtok_search = dimtok;
if ~isempty(cfg.features)
  dimtok_search(cfg.features-1) = [];
end

%% transform neighbours into boolean matrix if necessary
if isempty(cfg.mvpa.neighbours) && has_dim && has_dimord
  if ~isempty(cfg.timwin) && ~any(contains(dimtok_search, 'time'))
    ft_warning('ignoring cfg.timwin because time is not a search dimension')
  end
  if ~isempty(cfg.freqwin) && ~any(contains(dimtok_search, 'f'))
    ft_warning('ignoring cfg.freqwin because freq is not a search dimension')
  end
  if (~isempty(cfg.neighbours) || ~isempty(cfg.connectivity)) && ~any(contains(dimtok_search, 'chan'))
    ft_warning('ignoring cfg.connectivity and cfg.neighbours because chan is not a search dimension')
  end
  
  if (any(contains(dimtok_search, 'chan')) && (~isempty(cfg.neighbours) || ~isempty(cfg.connectivity))) || ...
      (any(contains(dimtok_search, 'time')) && ~isempty(cfg.timwin)) || ...
      (any(contains(dimtok_search, 'freq')) && ~isempty(cfg.freqwin))
    cfg.mvpa.neighbours = cell(numel(dimtok_search),1);
    for ix = 1:numel(dimtok_search)
      
      switch(dimtok_search{ix})
        case 'chan'
          % create boolean neighbour matrix for chan
          if isempty(cfg.neighbours) && ~isempty(cfg.connectivity)
            cfg.neighbours = cfg.connectivity;
          end
          if isstruct(cfg.neighbours)
            tmp_cfg = cfg;
            tmp_cfg.neighbours = cfg.neighbours;
            cfg.neighbours = channelconnectivity(tmp_cfg);
            cfg.neighbours = logical(double(cfg.neighbours) + eye(size(cfg.neighbours))); % include source channel
          end
          cfg.mvpa.neighbours{ix} = cfg.neighbours;
        case 'time'
          timdim = strcmp(dimtok, 'time');
          if ~isempty(cfg.timwin)
            % create boolean neighbour matrix for time
            T = ones(cfg.dim(timdim));
            cfg.mvpa.neighbours{ix} = T - triu(T, floor(cfg.timwin./2)+1) - tril(T, -floor(cfg.timwin./2)-1) > 0;
          else
            cfg.mvpa.neighbours{ix} = eye(cfg.dim(timdim));
          end
        case 'freq'
          % create boolean neighbour matrix for freq
          freqdim = strcmp(dimtok, 'freq');
          if ~isempty(cfg.freqwin)
            F = ones(cfg.dim(freqdim));
            cfg.mvpa.neighbours{ix} = F - triu(F, floor(cfg.freqwin./2)+1) - tril(F, -floor(cfg.freqwin./2)-1) > 0;
          else
            cfg.mvpa.neighbours{ix} = eye(cfg.dim(freqdim));
          end
        otherwise
          ft_error('Search dimension is ''%s'' but only ''chan'' ''freq'' and ''time'' are supported', dimtok_search{ix})
      end
    end
  end
elseif ~isempty(cfg.neighbours) || ~isempty(cfg.connectivity) || ~isempty(cfg.timwin) || ~isempty(cfg.freqwin)
  ft_warning('cfg.mvpa.neighbours has been set, ignoring cfg.neighbours/cfg.connectivity/cfg.timwin/cfg.freqwin')
end

%% adapt channel labels
if any(strcmp('chan', cfg.mvpa.dimension_names(cfg.features)))
  % combine all labels when chan is used as features
  label = sprintf('combined(%s)', strjoin(cfg.channel, ','));
elseif ~isempty(cfg.neighbours)
  label = cell(size(cfg.neighbours,1), 1);
  if (size(cfg.neighbours,1) == size(cfg.neighbours,2)) && all(diag(cfg.neighbours))
    % keep label unchanged
    label = cfg.channel;
  else
    selchan = find(strcmp(dimtok_search, 'chan'));
    % merge neighbours into combined channels
    for ix = 1:numel(label)
      chan_ix = cfg.mvpa.neighbours{selchan}(ix,:)>0;
      if sum(chan_ix)>1
        label{ix} = sprintf('combined(%s)', strjoin(cfg.channel(chan_ix), ','));
      else
        label{ix} = cfg.channel{chan_ix};
      end
    end
  end
end

%% Call MVPA-Light
if isempty(cfg.mvpa.model)
  % -------- Classification --------
  if ndims(dat)==3 && numel(cfg.features)==1 && cfg.features==2 && numel(cfg.generalize)==1 && cfg.generalize==3 && isempty(cfg.mvpa.neighbours)
    % special case: time generalization for 3D data
    [perf, result] = mv_classify_timextime(cfg.mvpa, dat, design);
  else
    [perf, result] = mv_classify(cfg.mvpa, dat, design);
  end
else
  % -------- Regression --------
  [perf, result] = mv_regress(cfg.mvpa, dat, design);
end

% build dimord from result struct
if has_dimord
  if ~iscell(perf)
    dimord = strrep(result.perf_dimension_names, ' ', '');
  else 
    % more than one output metric is requested, use the first one for the dimord
    dimord = strrep(result.perf_dimension_names{1}, ' ', '');
  end
  if iscell(dimord), dimord = strjoin(dimord, '_'); end
end

if ~iscell(cfg.mvpa.metric), cfg.mvpa.metric = {cfg.mvpa.metric}; end
if ~iscell(perf),            perf            = {perf};            end

%% setup stat struct
stat = [];
for mm=1:numel(perf)
  
  % Performance metric
  stat.(cfg.mvpa.metric{mm}) = perf{mm};
  
  % Std of performance
  if iscell(result.perf_std)
    stat.([cfg.mvpa.metric{mm} '_std']) = result.perf_std{mm};
  else
    stat.([cfg.mvpa.metric{mm} '_std']) = result.perf_std;
  end
end

% return the MVPA-Light result struct as well
stat.mvpa = result;

if isfield(cfg, 'latency') && ((isfield(cfg,'avgovertime') && strcmp(cfg.avgovertime, 'yes')) || (~isempty(cfg.mvpa.dimension_names) && any(ismember('time', cfg.mvpa.dimension_names(cfg.features)))))
  time = mean(cfg.latency);
end
if isfield(cfg, 'frequency')
  frequency = mean(cfg.frequency);
end

if exist('label', 'var'),     stat.label  = label;  end
if exist('dim', 'var'),       stat.dim    = dim;    end
if exist('dimord', 'var'),    cfg.dimord  = dimord; end % stat.dimord is overwritten by cfg.dimord in the caller, hence it's useless to set stat.dimord here
if exist('frequency', 'var'), stat.freq   = frequency; end
if exist('time', 'var'),      stat.time   = time; end

% helper functions
  function cfg = add_to_preprocess(cfg, item)
    if ~isfield(cfg.mvpa,'preprocess')
      cfg.mvpa.preprocess = item;
    elseif ~iscell(cfg.mvpa.preprocess)
      cfg.mvpa.preprocess = {item cfg.mvpa.preprocess};
    else
      cfg.mvpa.preprocess = [{item} cfg.mvpa.preprocess];
    end
  end
end
