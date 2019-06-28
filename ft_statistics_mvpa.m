function [stat, cfg] = ft_statistics_mvpa(cfg, dat, design)

% FT_STATISTICS_MVPA performs multivariate pattern classification
% on the data. If the data has not been averaged over time, classification
% is performed separately for every time point. Additionally, searchlight
% analysis (classification for each channel/voxel
% separately) or time x time generalisation can be performed.
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
%   cfg.searchlight     = 'yes' or 'no', performs searchlight analysis
%                         (default 'no'). More information see below
%   cfg.timextime       = 'yes' or 'no', performs time x time
%                         generalisation. In other words, the classifier
%                         is trained at each time point and tested at
%                         every time point. The result is a time x time
%                         matrix of classification performance.
%                         (default 'no')
%                         Note that searchlight and timextime cannot be
%                         run simultaneously (at least one option needs
%                         to be set to 'no').
%   cfg.mvpa            = structure that contains detailed options for the
%                         MVPA procedure. See
%                         https://github.com/treder/MVPA-Light for more
%                         details.
%
%
%   cfg.mvpa.classifier      = 'lda'          Regularised linear discriminant analysis
%                                        (LDA) (for two classes)
%                         'multiclass_lda' LDA for more than two classes
%                         'logreg'       Logistic regression
%                         'svm'          Support Vector Machine (SVM)
%                         'ensemble'     Ensemble of classifiers. Any of the other
%                                        classifiers can be used as a learner.
%                         'kernel_fda'   Kernel Fisher Discriminant Analysis
%   cfg.mvpa.metric          = string, performance metric. Possible metrics:
%                         accuracy auc tval dval confusion precision
%                         recall f1
%   See https://github.com/treder/MVPA-Light for an overview of all
%   classifiers and metrics.
%
%   cfg.mvpa.param           = struct, structure with hyperparameters for the
%                         classifier (see HYPERPARAMETERS below)
%
%   cfg.mvpa.balance         = string, for imbalanced data that does not have
%                         the same number of instances in each class
%                         'oversample'     oversamples the minority classes
%                         'undersample'    undersamples the minority classes
%                         such that all classes have the same number of
%                         instances. Note that undersampling is at the
%                         level of the repeats, whereas we oversampling occurs within each
%                         training set (for an explanation see mv_balance_classes).
%                         You can also give an integer number for undersampling.
%                         The samples will be reduced to this number. Note that
%                         concurrent over/undersampling (oversampling of the
%                         smaller class, undersampling of the larger class) is not
%                         supported at the moment
%  cfg.mvpa.replace          = bool, if balance is set to 'oversample' or 'undersample',
%                         replace determines whether data is drawn with
%                         replacement (default 1)
%  cfg.mvpa.normalise        = string, normalises the data across samples, for each time point
%                         and each feature separately, using 'zscore' or 'demean'
%                         (default 'zscore'). Set to 'none' or [] to avoid normalisation.
%  cfg.mvpa.feedback         = 'yes' or 'no', whether or not to print feedback on the console (default 'yes')
%
% To obtain a realistic estimate of classification performance,
% cross-validation is used. It is controlled by the following parameters:
%   cfg.mvpa.cv              = string, cross-validation type, either 'kfold', 'leaveout'
%                         or 'holdout'. If 'none', no cross-validation is
%                         used and the classifier is tested on the training
%                         set. (default 'kfold')
%   cfg.mvpa.k               = number of folds in k-fold cross-validation (default 5)
%   cfg.mvpa.repeat          = number of times the cross-validation is repeated
%                         with new randomly assigned folds (default 5)
%   cfg.mvpa.p               = if cfg.cv is 'holdout', p is the fraction of test
%                         samples (default 0.1)
%   cfg.mvpa.stratify        = if 1, the class proportions are approximately
%                         preserved in each test fold (default 1)
%
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
% Classification using a feature searchlight approach can highlight which
% feature(s) are informative. To this end, classification is performed on
% each feature separately. However, neighbouring features can enter the
% classification together when specified. Optional parameters:
%
% cfg.mvpa.neighbours   = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%                         Alternatively, a [features x features] matrix specifying
%                         which features are neighbours of each other. This
%                         matrix can be either
%                         a GRAPH consisting of 0's and 1's. A 1 in the
%                         (i,j)-th element signifies that feature i and feature j
%                         are neighbours, and a 0 means they are not neighbours
%                         or a DISTANCE MATRIX, where larger values mean larger distance.
%                         If no matrix is provided, every feature is only neighbour
%                         to itself and classification is performed for each feature
%                         separately.
% cfg.mvpa.size         = size of the 'neighbourhood' of a feature.
%                         number of steps taken through the neighbourhood
%                         to include neighbours
%                         0: only the feature itself is considered (no neighbours)
%                         1: the feature and its immediate neighbours
%                         2: the feature, its neighbours, and its neighbours'
%                         neighbours
%                         3+: neighbours of neighbours of neighbours etc
%                         (default 1)
%                         if cfg.neighbours is a distance matrix, size defines the number of
%                         neighbouring features that enter the classification
%                         0: only the feature itself is considered (no neighbours)
%                         1: the feature and its first closest neighbour
%                            according to the distance matrix
%                         2+: the 2 closest neighbours etc.
%
% -- TODO: for time x time generalisation, in MVPA light we can use two
% different datasets (one for training the classifier, the other one for
% testing). This could be realised eg using extra fields in cfg such as
% cfg.X2 for dataset 2 and cfg.design2 for the second design matrix
%
% Returns:
%   stat        = struct with results. the .metric field contains the
%                   requested metrics
%

% Copyright (C) 2019, Matthias Treder
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

% cfg: set defaults
cfg.searchlight     = ft_getopt(cfg, 'searchlight', 'no');
cfg.timextime       = ft_getopt(cfg, 'timextime',   'no');
cfg.mvpa            = ft_getopt(cfg, 'mvpa',        []);
cfg.mvpa.metric     = ft_getopt(cfg.mvpa, 'metric',   'accuracy');
cfg.mvpa.feedback   = ft_getopt(cfg.mvpa, 'feedback',   'yes');
cfg.mvpa.neighbours = ft_getopt(cfg.mvpa, 'neighbours', []);

% flip dimensions such that the number of trials comes first
dat = dat';

% if cfg.dim has two entries which are non-singleton then the data has been
% 3D before. We must reshape it from 2D to 3D to run it in MVPA-Light
data_is_3D = (numel(cfg.dim) > 1 && all(cfg.dim > 1));  % checks whether data was 3D before being reshaped
if data_is_3D
  dat = reshape(dat, size(dat,1), cfg.dim(1), cfg.dim(2));
end

y = design;

%% perform sanity checks on parameters
if istrue(cfg.timextime) && istrue(cfg.searchlight)
  ft_error('you should not set timextime = ''yes'' and searchlight = ''yes'' simultaneously')
end

% timextime = 'yes' but data is not 3D we should change timextime to 'no'
if istrue(cfg.timextime) && ~data_is_3D
  ft_warning('timextime = ''yes'' but data has no time dimension, setting timextime = ''no''');
  cfg.timextime = 'no';
end

label = [];
dim   = [];
dimord = [];

%% Call MVPA-Light

if strcmp(cfg.searchlight, 'yes')
  % --- searchlight analysis ---

  if isstruct(cfg.mvpa.neighbours)
    cfg.mvpa.neighbours = channelconnectivity(struct('neighbours',cfg.mvpa.neighbours, 'channel', {cfg.channel}));
  end
  [perf, result] = mv_searchlight(cfg.mvpa, dat, y);
  
  % this preserves any spatial dimension, so no adjustment is done to a
  % channel list, if present
  if isfield(cfg, 'channel')
    label = cfg.channel;
  end

  if isfield(cfg, 'dim')
    dim = cfg.dim;
  end

elseif strcmp(cfg.timextime, 'yes')
  % --- time x time generalisation ---
  [perf, result] = mv_classify_timextime(cfg.mvpa, dat, y);

  % this does note preserve any spatial dimension, so label should be
  % adjusted
  label = squeezelabel(label, cfg);
  dim   = squeezedim(dim, cfg);
  dimord = 'time_time';

elseif data_is_3D
  % --- classification across time ---
  [perf, result] = mv_classify_across_time(cfg.mvpa, dat, y);

  % this does note preserve any spatial dimension, so label should be
  % adjusted
  label = squeezelabel(label, cfg);
  dim   = squeezedim(dim, cfg);

else
  % --- data has no time dimension, perform only cross-validation ---
  [perf, result] = mv_crossvalidate(cfg.mvpa, dat, y);

  % this does note preserve any spatial dimension, so label should be
  % adjusted
  label = squeezelabel(label, cfg);
  dim   = squeezedim(dim, cfg);

end

%% setup stat struct
stat = [];
if ~iscell(cfg.mvpa.metric), cfg.mvpa.metric = {cfg.mvpa.metric}; end
if ~iscell(perf),            perf            = {perf};            end
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

if ~isempty(label)
  stat.label = label;
end

if ~isempty(dim)
  stat.dim = dim;
end

if ~isempty(dimord)
  stat.dimord = dimord;
end

function label = squeezelabel(label, cfg)

if isfield(cfg, 'channel')
  label = sprintf('combined(%s)', sprintf('%s',cfg.channel{:}));
end

function dim = squeezedim(dim, cfg)

if isfield(cfg, 'dim')
  dim = cfg.dim;
  dim(1) = 1;
end
