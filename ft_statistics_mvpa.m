function [stat, cfg] = ft_statistics_mvpa(cfg, dat, design)

% FT_STATISTICS_MVPA performs multivariate pattern classification or 
% regression using the MVPA-Light toolbox. The function supports
% cross-validation, searchlight analysis, generalization, nested
% preprocessing, a variety of classification and regression metrics, as
% well as statistical testing of these metrics.
%
% This function should not be called directly, instead
% you should call the function that is associated with the type of
% data on which you want to perform the test.
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
% where the data is obtained from FT_TIMELOCKANALYSIS, FT_FREQANALYSIS
% or FT_SOURCEANALYSIS respectively, or from FT_TIMELOCKGRANDAVERAGE,
% FT_FREQGRANDAVERAGE or FT_SOURCEGRANDAVERAGE respectively.
%
% The configuration can contain
%   cfg.features        = specifies the name or index of the dimension(s) 
%                         that serve(s) as features for the classifier or
%                         regression model. Dimensions that are not
%                         samples or features act as searchlight
%                         dimensions. For instance, assume the data is a
%                         3D array of size [samples x channels x time]. 
%                         If cfg.features = 2, the channels serve as 
%                         features. A classification is then performed for
%                         each time point (we call time a searchlight
%                         dimension). Conversely, if cfg.features = 3, the
%                         time points serve as features. A classification
%                         is performed for each channel (channel is a
%                         searchlight dimension). 
%                         If cfg.features = [], then all non-sample
%                         dimensions serve as searchlight dimensions.
%                         If the dimensions have names (ie cfg.dimord
%                         exists), then instead of numbers the feature can
%                         be specified as a string (e.g. 'chan').
%                         (default 2)
%   cfg.generalize      = specifies the name or index of the dimensions
%                         that serves for generalization (if any). For
%                         instance, if the data is [samples x channels x
%                         time], and cfg.generalize = 3, a time x time
%                         generalization is performed. If cfg.generalize =
%                         2, a electrode x electrode generalization is
%                         performed. cfg.generalize must refer to a
%                         searchlight dimension, therefore its value must
%                         be different from the value of cfg.features.
%                         (default [])
%   cfg.mvpa            = structure that contains detailed options for the
%                         MVPA procedure. See
%                         https://github.com/treder/MVPA-Light for more
%                         details on the parameters and the available
%                         statistical models and metrics.
%
%   cfg.mvpa.classifier = string specifying the classifier 
%                         Available classifiers:
%                         'ensemble'     Ensemble of classifiers. Any of the other
%                                        classifiers can be used as a learner.
%                         'kernel_fda'   Kernel Fisher Discriminant Analysis
%                         'lda'          Regularised linear discriminant analysis
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
%   cfg.mvpa.model      = string specifying the regression model. If a
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
%   cfg.mvpa.metric     = string, classification or regression metric. 
%                         Classification metrics: accuracy auc confusion 
%                         dval f1 kappa precision recall tval
%                         Regression metrics: mae mse r_squared
%                         
%   cfg.mvpa.hyperparameter  = struct, structure with hyperparameters for the
%                         classifier or regression model (see HYPERPARAMETERS below)
%   cfg.mvpa.feedback   = 'yes' or 'no', whether or not to print feedback on the console (default 'yes')
%
% To obtain a realistic estimate of classification performance,
% cross-validation is used. It is controlled by the following parameters:
%   cfg.mvpa.cv         = string, cross-validation type, either 'kfold', 'leaveout'
%                         'holdout', or 'predefined'. If 'none', no cross-validation is
%                         used and the model is tested on the training
%                         set. (default 'kfold')
%   cfg.mvpa.k          = number of folds in k-fold cross-validation (default 5)
%   cfg.mvpa.repeat     = number of times the cross-validation is repeated
%                         with new randomly assigned folds (default 5)
%   cfg.mvpa.p          = if cfg.cv is 'holdout', p is the fraction of test
%                         samples (default 0.1)
%   cfg.mvpa.stratify   = if 1, the class proportions are approximately
%                         preserved in each test fold (default 1)
%   cfg.mvpa.fold       = if cv='predefined', fold is a vector of length
%                         #samples that specifies the fold each sample belongs to
%
% More information about each classifier is found in the documentation of
% MVPA-Light (github.com/treder/MVPA-Light/).
%
% HYPERPARAMETERS:
% Each classifier comes with its own set of hyperparameters, such as
% regularisation parameters and the kernel. Hyperparameters can be set
% using the cfg.param substruct. For instance, in LDA, cfg.param.lambda =
% 'auto' sets the lambda regularisation parameter.
%
% The specification of the hyperparameters is found in the training function
% for each classifier at github.com/treder/MVPA-Light/tree/master/classifier
% If a hyperparameter is not specified, default values are used.
%
% SEARCHLIGHT ANALYSIS:
% Data dimensions that are not samples or features serve as 'search
% dimensions'. For instance, if the data is [samples x chan x time] 
% and cfg.features = 'time', then the channel dimension serves as search
% dimension: a separate analysis is carried out for each channel. Instead
% of considering each channel individually, a searchlight can be defined 
% such that each channel is used together with its neighbours. Neighbours
% can be specified using the cfg.mvpa.neighbours field:
%
%   cfg.mvpa.neighbours   = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%                           Alternatively, a [features x features] matrix specifying
%                           which features are neighbours of each other. This
%                           matrix consists of 0's and 1's. A 1 in the
%                           (i,j)-th element signifies that feature i and feature j
%                           are neighbours, and a 0 means they are not neighbours.
%
% TODO: for time x time generalisation, in MVPA light we can use two
% different datasets (one for training the classifier, the other one for
% testing). This could be realised eg using extra fields in cfg such as
% cfg.X2 for dataset 2 and cfg.design2 for the second design matrix
%
% Returns:
%   stat        = struct with results. the .metric field contains the
%                   requested metrics
%

% Copyright (C) 2019-2020, Matthias Treder and Jan-Mathijs Schoffelen
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
cfg.features        = ft_getopt(cfg, 'features',    2);
cfg.generalize      = ft_getopt(cfg, 'generalize',  []);
cfg.connectivity    = ft_getopt(cfg, 'connectivity', []); % the default is dealt with below

cfg.mvpa            = ft_getopt(cfg, 'mvpa',         []);
cfg.mvpa.neighbours = ft_getopt(cfg.mvpa, 'neighbours', []);
cfg.mvpa.model      = ft_getopt(cfg.mvpa, 'model', []);
if isempty(cfg.mvpa.model)
  cfg.mvpa.classifier = ft_getopt(cfg.mvpa, 'classifier', 'lda');
  if strcmp(cfg.mvpa.classifier, 'naive_bayes')
    cfg.mvpa.append = ft_getopt(cfg.mvpa, 'append', 1);
  end
  cfg.mvpa.metric     = ft_getopt(cfg.mvpa, 'metric', 'accuracy');
else
  cfg.mvpa.metric     = ft_getopt(cfg.mvpa, 'metric', 'mae');
end
cfg.mvpa.feedback   = ft_getopt(cfg.mvpa, 'feedback',   'yes');

if isfield(cfg, 'dimord')
  cfg.mvpa.dimension_names = ft_getopt(cfg.mvpa, 'dimension_names', [{'samples'} tokenize(cfg.dimord, '_')]);
end

% flip dimensions such that the number of trials comes first
dat = dat.';

if isfield(cfg, 'dim')
  % MVPA-Light expects the original multi-dimensional array
  dat = reshape(dat, [size(dat,1) cfg.dim]);
end

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

%% convert features and generalize from char to dimension numbers
if ischar(cfg.features)
  assert(isfield(cfg, 'dimord'), 'if cfg.features is a string then cfg.dimord must exist')
  cfg.features = find(ismember(cfg.mvpa.dimension_names, cfg.features));
  if isempty(cfg.features)
    ft_error(sprintf('cfg.features = ''%s'' is not contained in cfg.dimord', cfg.features))
  end
end

if ischar(cfg.generalize)
  assert(isfield(cfg, 'dimord'), 'if cfg.generalize is a string then cfg.dimord must exist')
  cfg.generalize = find(ismember(cfg.mvpa.dimension_names, cfg.generalize));
  if isempty(cfg.generalize)
    ft_error(sprintf('cfg.generalize = ''%s'' is not contained in cfg.dimord', cfg.generalize))
  end
end

cfg.mvpa.feature_dimension          = cfg.features;
cfg.mvpa.generalization_dimension   = cfg.generalize;

if any(strcmp('chan', cfg.mvpa.dimension_names(cfg.features)))
  label = sprintf('combined(%s)', sprintf('%s',cfg.channel{:})); % combine labels when chan is used as features
end

%% transform neighbours into connectivity matrix
if ~isfield(cfg.mvpa,'neighbours')
  if ~isempty(cfg.connectivity)
    cfg.mvpa.neighbours = cfg.connectivity;
  else
    cfg.mvpa.neighbours = [];
  end
end

if isstruct(cfg.mvpa.neighbours)
  tmp_cfg = cfg;
  tmp_cfg.neighbours = cfg.mvpa.neighbours;
  cfg.mvpa.neighbours = channelconnectivity(tmp_cfg);
elseif iscell(cfg.mvpa.neighbours)
  tmp_cfg = cfg;
  for ix = 1:numel(cfg.mvpa.neighbours)
    tmp_cfg.neighbours{ix} = cfg.mvpa.neighbours{ix};
    cfg.mvpa.neighbours{ix} = channelconnectivity(tmp_cfg);
  end
end

%% Call MVPA-Light
if isempty(cfg.mvpa.model)
  % -------- Classification --------
  if ndims(dat)==3 && cfg.mvpa.feature_dimension==2 && cfg.mvpa.generalization_dimension==3
    [perf, result] = mv_classify_timextime(cfg.mvpa, dat, design);
  else
    [perf, result] = mv_classify(cfg.mvpa, dat, design);
  end
else
  % -------- Regression --------
  [perf, result] = mv_regress(cfg.mvpa, dat, design);
end

% build dimord from result struct
if isfield(cfg, 'dimord')
  dimord = strrep(result.perf_dimension_names, ' ', '');
  if iscell(dimord), dimord = strjoin(dimord, '_'); end
end

if ~iscell(cfg.mvpa.metric), cfg.mvpa.metric = {cfg.mvpa.metric}; end
if ~iscell(perf),            perf            = {perf};            end

%% check which data dim descriptors need to be updated
if isfield(cfg, 'latency')
  time = mean(cfg.latency);
end
if isfield(cfg, 'frequency')
  frequency = mean(cfg.frequency);
end

%% setup stat struct
stat = [];

% return the MVPA-Light result struct as well
stat.mvpa = result;

if exist('label', 'var'),     stat.label  = label;  end
if exist('dim', 'var'),       stat.dim    = dim;    end
if exist('dimord', 'var'),    stat.dimord = dimord; end
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
