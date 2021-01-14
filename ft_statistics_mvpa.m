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
%   cfg.feature         = specifies the name or index of the dimension 
%                         that serves as features for the classifier or
%                         regression model. Dimensions that are not
%                         samples or features act as searchlight
%                         dimensions. For instance, assume the data is a
%                         3D array of size [samples x channels x time]. 
%                         If cfg.feature = 2, the channels serve as 
%                         features. A classification is then performed for
%                         each time point (we call time a searchlight
%                         dimension). Conversely, if cfg.feature = 3, the
%                         time points serve as features. A classification
%                         is performed for each channel (channel is a
%                         searchlight dimension). 
%                         If cfg.feature = [], then all non-sample
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
%                         be different from the value of cfg.feature.
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
%
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
%                         used and the classifier is tested on the training
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
% HYPERPARAMETERS:
% Each classifier comes with its own set of hyperparameters, such as
% regularization parameters and the kernel. Hyperparameters can be set
% using the cfg.hyperparameter substruct. For instance, in LDA, 
% cfg.hyperparameter.lambda = 'auto' sets the lambda regularization parameter.
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

% backward compatibility
cfg.mvpa = ft_checkconfig(cfg.mvpa, 'renamed', {'param', 'hyperparameter'});
ft_checkconfig(cfg.mvpa, 'deprecated', {'balance' 'normalise' 'replace'});

% cfg: set defaults 
cfg.feature         = ft_getopt(cfg, 'feature',     2);
cfg.generalize      = ft_getopt(cfg, 'generalize',  []);
cfg.connectivity    = ft_getopt(cfg, 'connectivity', []); % the default is dealt with below

cfg.mvpa            = ft_getopt(cfg, 'mvpa',         []);
cfg.mvpa.neighbours = ft_getopt(cfg, 'neighbours', []);
cfg.mvpa.classifier = ft_getopt(cfg.mvpa, 'classifier', []);
cfg.mvpa.model      = ft_getopt(cfg.mvpa, 'model', []);
cfg.mvpa.metric     = ft_getopt(cfg.mvpa, 'metric',     'accuracy');
cfg.mvpa.feedback   = ft_getopt(cfg.mvpa, 'feedback',   'yes');
cfg.mvpa.timwin     = ft_getopt(cfg.mvpa, 'timwin',     []); % in samples because time axis might not be known

assert(isempty(cfg.mvpa.classifier) || isempty(cfg.mvpa.model), 'set either a cfg.mvpa.classifier (for classification) or a cfg.mvpa.model (for regression), not both')

if isfield(cfg, 'dimord')
    cfg.mvpa.dimension_names = ft_getopt(cfg.mvpa, 'dimension_names', [{'samples'} tokenize(cfg.dimord, '_')]);
end

% deal with the neighbourhood of the channels/triangulation/voxels
if isempty(cfg.connectivity)
  if isfield(cfg, 'dim') && ~isfield(cfg, 'channel') && ~isfield(cfg, 'tri')
    % input data can be reshaped into a 3D volume, use bwlabeln/spm_bwlabel rather than clusterstat
    ft_info('using connectivity of voxels in 3-D volume\n');
    cfg.connectivity = nan;
    %if isfield(cfg, 'inside')
    %  cfg = fixinside(cfg, 'index');
    %end
  elseif isfield(cfg, 'tri')
    % input data describes a surface along which neighbours can be defined
    ft_info('using connectivity of vertices along triangulated surface\n');
    cfg.connectivity = triangle2connectivity(cfg.tri);
    if isfield(cfg, 'insideorig')
      cfg.connectivity = cfg.connectivity(cfg.insideorig, cfg.insideorig);
    end
  elseif isfield(cfg, 'avgoverchan') && istrue(cfg.avgoverchan)
    % channel dimension has been averaged across, no sense in clustering across space
    cfg.connectivity = true(1);
  elseif isfield(cfg, 'channel')
    cfg.neighbours   = ft_getopt(cfg, 'neighbours', []);
    cfg.connectivity = channelconnectivity(cfg);
  else
    % there is no connectivity in the spatial dimension
    cfg.connectivity = false(size(dat,1));
  end
else
  % use the specified connectivity: op hoop van zegen
end

if isfield(cfg, 'dim') && isfield(cfg, 'dimord') && contains(cfg.dimord, 'time') && ~isempty(cfg.mvpa.timwin)
  % create neighourhood matrix for time dimension
  dimtok = tokenize(cfg.dimord, '_');
  timdim = find(strcmp(dimtok, 'time'));
  T = ones(cfg.dim(timdim));
  T = T - triu(T, floor(cfg.mvpa.timwin./2)) - tril(T, -ceil(cfg.mvpa.timwin./2)) > 0;
  if isfield(cfg.mvpa, 'time')&&~isempty(cfg.mvpa.time)
    % this cannot be dealt with by ft_getopt, because the lower level code
    % will in that case not properly set the default value
    T = T(cfg.mvpa.time,:);
  end
else
  T = [];
end

% flip dimensions such that the number of trials comes first
dat = dat.';

% MVPA-Light expects the original multi-dimensional array
dat = reshape(dat, [size(dat,1) cfg.dim]);

%% convert features and generalize from char to dimension numbers
if ischar(cfg.feature)
    assert(isfield(cfg, 'dimord'), 'if cfg.feature is a string then cfg.dimord must exist')
    cfg.feature = find(ismember(cfg.mvpa.dimension_names, cfg.feature));
    if isempty(cfg.feature)
        ft_error(sprintf('cfg.feature = ''%s'' is not found in cfg.dimord', cfg.feature))
    end
end

if any(strcmp('chan',cfg.mvpa.dimension_names(cfg.feature)))
    % combine labels
    label = sprintf('combined(%s)', sprintf('%s',cfg.channel{:}));
end

if ischar(cfg.generalize)
    assert(isfield(cfg, 'dimord'), 'if cfg.generalize is a string then cfg.dimord must exist')
    cfg.generalize = find(ismember(cfg.mvpa.dimension_names, cfg.generalize));
    if isempty(cfg.generalize)
        ft_error(sprintf('cfg.generalize = ''%s'' is not contained in cfg.dimord', cfg.generalize))
    end
end

cfg.mvpa.feature_dimension          = cfg.feature;
cfg.mvpa.generalization_dimension   = cfg.generalize;

%% Call MVPA-Light
if isempty(cfg.mvpa.model)
    % -------- Classification --------
    [perf, result] = mv_classify(cfg.mvpa, dat, design);
else
    % -------- Regression --------
    [perf, result] = mv_regress(cfg.mvpa, dat, design);
end

% build dimord from result struct
dimord = strrep(result.perf_dimension_names, ' ', '');

% if istrue(cfg.searchlight) && isempty(T)
%   % --- searchlight analysis across the spatial dimension ---
% 
%   if ~isequal(cfg.connectivity, false(size(cfg.connectivity))) && ~isempty(cfg.connectivity)
%     cfg.mvpa.neighbours = cfg.connectivity;
%   end
%   [perf, result] = mv_searchlight(cfg.mvpa, dat, y);
%   
%   % create boolean vector for the update of dimension descriptors at higher
%   % level structure
%   adj_dim = cfg.dim ~= size(perf); % assumes first element to be spatial dimension
%   
% elseif istrue(cfg.searchlight) && ~isempty(T) && data_is_3D
%   % --- searchlight across time, or searchlight across space and time ---
%   
%   if ~isequal(cfg.connectivity, false(size(cfg.connectivity))) && ~isempty(cfg.connectivity)
%     cfg.mvpa.neighbours        = {cfg.connectivity T};
%     cfg.mvpa.sample_dimension  = 1;
%     cfg.mvpa.feature_dimension = [];
%     if isequal(cfg.mvpa.classifier, 'naive_bayes')
%       cfg.mvpa.append = true;
%     end
%   else
%     cfg.mvpa.neighbours        = T;
%     cfg.mvpa.sample_dimension  = 1;
%     cfg.mvpa.feature_dimension = 2; 
%   end
%   
%   [perf, result] = mv_classify(cfg.mvpa, dat, y);
%   
%   if iscell(cfg.mvpa.neighbours)
%     adj_dim = cfg.dim ~= size(perf);
%   else 
%     adj_dim = [true false(1,numel(cfg.dim)-1)];
%   end
%   
%   if numel(perf)~=prod(cfg.dim) && all(~adj_dim)
%     dimorig = cfg.dim;
%     cfg.dim = [size(cfg.mvpa.hyperparameter.neighbours{1},1), ...
%                size(cfg.mvpa.hyperparameter.neighbours{2},1)];
%     perf    = reshape(perf, cfg.dim);
%     result.perf_std = reshape(result.perf_std, cfg.dim);
%     
%     if dimorig(1)~=cfg.dim(1)
%       label = cell(cfg.dim(1),1);
%       for k = 1:cfg.dim(1)
%         label{k} = sprintf('feature%04d',k);
%       end
%     end
%     if dimorig(2)~=cfg.dim(2)&&cfg.dim(2)&&isfield(cfg,'latency')
%       timeorig = linspace(cfg.latency(1),cfg.latency(2),size(T,2));
%       time = zeros(1,size(T,1));
%       for k = 1:numel(time)
%         time(k) = mean(timeorig(T(k,:)));
%       end
%     end
%   end
%   
% elseif istrue(cfg.timextime)
%   % --- time x time generalisation ---
%   [perf, result] = mv_classify_timextime(cfg.mvpa, dat, y);
%   
%   % this does note preserve any spatial dimension, so label should be
%   % adjusted
%   label = sprintf('combined(%s)', sprintf('%s',cfg.channel{:}));
%   dim   = squeezedim(cfg.dim);
%   dimord = 'time_time';
%   adj_dim = [false false];
%   
% elseif data_is_3D && ~istrue(cfg.searchlight)
%   % --- classification across the non-spatial dimension ---
%   [perf, result] = mv_classify_across_time(cfg.mvpa, dat, y);
%   
%   adj_dim = [true false(1,numel(cfg.dim)-1)];
%   
% else
%   % this is the generic fallback, which seems very similar to the second
%   % instance above...
%   if ~isempty(T) && ~isequal(cfg.connectivity, false(size(cfg.connectivity))) && ~isempty(cfg.connectivity)
%     cfg.mvpa.hyperparameter.neighbours = {cfg.connectivity T};
%     adj_dim = false(size(cfg.dim));
%   elseif ~isempty(T)
%     cfg.mvpa.hyperparameter.neighbours = {eye(cfg.dim(1)) T};
%     adj_dim = false(size(cfg.dim));
%   else
%     adj_dim = true(size(cfg.dim));
%   end
%   [perf, result] = mv_crossvalidate(cfg.mvpa, dat, y);
%   
%   if numel(perf)~=prod(cfg.dim) && all(~adj_dim)
%     dimorig = cfg.dim;
%     cfg.dim = [size(cfg.mvpa.hyperparameter.neighbours{1},1), ...
%                size(cfg.mvpa.hyperparameter.neighbours{2},1)];
%     perf    = reshape(perf, cfg.dim);
%     result.perf_std = reshape(result.perf_std, cfg.dim);
%     
%     if dimorig(1)~=cfg.dim(1)
%       label = cell(cfg.dim(1),1);
%       for k = 1:cfg.dim(1)
%         label{k} = sprintf('feature%04d',k);
%       end
%     end
%     if dimorig(2)~=cfg.dim(2)&&cfg.dim(2)&&isfield(cfg,'latency')
%       timeorig = linspace(cfg.latency(1),cfg.latency(2),size(T,2));
%       time = zeros(1,size(T,1));
%       for k = 1:numel(time)
%         time(k) = mean(timeorig(T(k,:)));
%       end
%     end
%   end
% end


if ~iscell(cfg.mvpa.metric), cfg.mvpa.metric = {cfg.mvpa.metric}; end
if ~iscell(perf),            perf            = {perf};            end

%% check which data dim descriptors need to be updated 
shiftflag = zeros(1,numel(perf));
if ~exist('label', 'var') && (isfield(cfg, 'channel') && adj_dim(1))
  label = sprintf('combined(%s)', sprintf('%s',cfg.channel{:}));
  % for consistency higher up, the first dimension of perf should be
  % singleton
  for k = 1:numel(perf)
    siz = size(perf{k});
    if siz(end) == 1 && siz(1) ~=1
      shiftflag(k) = -1;
    end
  end
end

if ~exist('dim', 'var') && isfield(cfg, 'dim') 
  dim = size(perf); % FIXME check whether this is always correct, I doubt it
end

if isfield(cfg, 'frequency') && ~adj_dim(2)
  if isfield(cfg, 'latency') && adj_dim(3)
    time = mean(cfg.latency);
  end
elseif isfield(cfg, 'frequency') && adj_dim(2)
  frequency = mean(cfg.frequency);
  if isfield(cfg, 'latency') && adj_dim(3)
    time = mean(cfg.latency);
  end
elseif ~isfield(cfg, 'frequency')
  if isfield(cfg, 'latency') && adj_dim(2)
    time = mean(cfg.latency);
  end
end

%% setup stat struct
stat = [];
for mm=1:numel(perf)

  % Performance metric
  stat.(cfg.mvpa.metric{mm}) = shiftdim(perf{mm}, shiftflag);

  % Std of performance
  if iscell(result.perf_std)
    stat.([cfg.mvpa.metric{mm} '_std']) = shiftdim(result.perf_std{mm}, shiftflag);
  else
    stat.([cfg.mvpa.metric{mm} '_std']) = shiftdim(result.perf_std, shiftflag);
  end
end

% return the MVPA-Light result struct as well
stat.mvpa = result;

if exist('label', 'var'),     stat.label  = label;  end
if exist('dim', 'var'),       stat.dim    = dim;    end
if exist('dimord', 'var'),    stat.dimord = dimord; end
if exist('frequency', 'var'), stat.freq   = frequency; end
if exist('time', 'var'),      stat.time   = time; end



function dim = squeezedim(dim)

dim(1) = 1;

