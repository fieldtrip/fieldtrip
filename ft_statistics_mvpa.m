function [stat, cfg] = ft_statistics_mvpa(cfg, dat, design)
% FT_STATISTICS_MVPA performs multivariate pattern classification
% on the data. A search parameter can be set to specify which data
% dimensions so search/loop over. Additionally, generalization (time x time
% or freq x freq) is implemented. 
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
%   cfg.search          = char or cell array, specifies which dimensions to
%                         search/loop over when the data is multi-dimensional. 
%                         If cfg.search = 'time' 
%                         then a separate classification is performed for
%                         every time point. If cfg.search = 'chan', a
%                         separate classification is performed for every
%                         channel. If cfg.search = 'freq', a separate
%                         classification is performed for every frequency.
%                         Multiple search dimensions can be given, e.g.
%                         {'time' 'freq'}. 
%                         Note that all data dimensions that are not search
%                         dimensions will be considered as features. E.g.,
%                         if the data is rpt x chan x time, and
%                         cfg.search='time', then the channels will be used
%                         as features, and if cfg.search='chan', the times
%                         will be used as features, unless either
%                         cfg.avgovertime or cfg.avgoverchan are true.
%                         In the search, neighbouring
%                         channels/times/frequencies can be added as features, 
%                         see below for more information on how to do this.
%                         (default values is 'time' if the data has a time 
%                         dimension, else 'freq' if the data has a frequency
%                         dimension, else '')
%   cfg.generalize      = 'no' or 'time' or 'freq'. If 'time', performs time
%                         generalization (aka time x time classification). 
%                         The classifier is trained at each time point and 
%                         tested at every time point. The result is a 
%                         time x time matrix of classification performance.
%                         For frequency or time-frequency data, it can be
%                         alternatively set to 'freq', so that
%                         generalization is performed across frequencies.
%                         (default 'no')
%   cfg.mvpa            = structure that contains detailed options for the
%                         MVPA procedure. See
%                         https://github.com/treder/MVPA-Light for more
%                         details.
%   cfg.mvpa.classifier      = 'lda'          Regularised linear discriminant analysis
%                                        (LDA) (for two classes)
%                         'multiclass_lda' LDA for more than two classes
%                         'logreg'       Logistic regression
%                         'naive_bayes'  Naive Bayes
%                         'svm'          Support Vector Machine (SVM)
%                         'ensemble'     Ensemble of classifiers. Any of the other
%                                        classifiers can be used as a learner.
%                         'kernel_fda'   Kernel Fisher Discriminant Analysis
%                          'libsvm'      (requires LIBSVM toolbox)
%                          'liblinear'   (requires LIBLINEAR toolbox)
%                          Documentation for each classifier is found at 
%                          github.com/treder/MVPA-Light/#classifiers
%   cfg.mvpa.metric          = string, performance metric. Possible metrics:
%                          accuracy auc tval dval confusion kappa precision
%                          recall f1. Multiple metrics can be given as a
%                          cell array eg {'accuracy' 'f1'}.
%                          Documentation for each metric is found at 
%                          github.com/treder/MVPA-Light/#metrics
%   cfg.mvpa.hyperparameter = struct, structure with hyperparameters for the
%                         classifier (see HYPERPARAMETERS below)
%   cfg.mvpa.feedback         = 'yes' or 'no', whether or not to print feedback on the console (default 'yes')
%
% CROSS-VALIDATION is used to obtain a realistic estimate of classification 
% performance. It is controlled by the following parameters:
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
% PREPROCESSING:  
% Defines a nested preprocessing pipeline. The preprocessing is applied
% on the traindata within the cross-validation. Parameters estimated
% from the train data (e.g. mean, variance) are then used to transform the
% test data. (see https://github.com/treder/MVPA-Light#preprocessing)
% .preprocess         - cell array containing the preprocessing pipeline. The
%                       pipeline is applied in chronological order
% .preprocess_param   - cell array of preprocessing parameter structs for each
%                       function. Length of preprocess_param must match length
%                       of preprocess
%
% HYPERPARAMETERS:
% Each classifier comes with its own set of hyperparameters, such as
% regularisation parameters and the kernel. Hyperparameters can be set
% using the cfg.hyperparameter substruct. For instance, in LDA, 
%         cfg.hyperparameter.lambda = 'auto' 
% sets the lambda regularisation parameter to automatic regularisation.
%
% The specification of the hyperparameters is found in the training function
% for each classifier at github.com/treder/MVPA-Light/tree/master/classifier
% If a hyperparameter is not specified, default values are used.
%
% SEARCH using NEIGHBOURS: 
% When one or more search dimensions are given, a separate
% classification is performed for each element along the search dimension
% (e.g. time, channel, or frequency). Such an analysis highlights which
% times/channels/frequencies contain discriminative information. 
% Sometimes it is useful to use a "searchlight" that includes neighbouring
% channels, time points, of frequencies, as well. To this end, the
% following optional parameters can be specified:
%
% cfg.mvpa.neighbours   = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%                         The neighbours must correspond to the search
%                         dimension. For instance, if cfg.search = 'chan',
%                         then the neighbours must correspond to channels.
%                         If cfg.search = 'time', they must correspond to
%                         the time points. If no neighbours are provided, 
%                         Alternatively, a matrix can be given whose size 
%                         matrix must agree with the search dimension (eg
%                         [channel x channel] or [time x time] matrix).
%                         This matrix can be either
%                         a GRAPH consisting of 0's and 1's. A 1 in the
%                         (i,j)-th element signifies that feature i and feature j
%                         are neighbours, and a 0 means they are not neighbours
%                         or a DISTANCE MATRIX, where larger values mean larger distance.
%                         If no matrix is provided, every feature is only neighbour
%                         to itself and classification is performed for each feature
%                         separately.
%                         Note that multiple neighbour structures can be
%                         given as a cell array if cfg.search contains
%                         multiple entries.
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
cfg.generalize      = ft_getopt(cfg, 'generalize',   'no');
cfg.mvpa            = ft_getopt(cfg, 'mvpa',        struct());
cfg.mvpa.classifier = ft_getopt(cfg.mvpa, 'classifier', 'lda');
cfg.mvpa.metric     = ft_getopt(cfg.mvpa, 'metric',     'accuracy');
cfg.mvpa.feedback   = ft_getopt(cfg.mvpa, 'feedback',   'yes');
cfg.mvpa.neighbours = ft_getopt(cfg.mvpa, 'neighbours', []);
cfg.mvpa.cv         = ft_getopt(cfg.mvpa, 'cv', 'kfold');

% check type of cfg parameters
ft_checkopt(cfg,      'generalize', 'char', {'no', 'time', 'freq'});
ft_checkopt(cfg,      'mvpa',       'struct');
ft_checkopt(cfg.mvpa, 'feedback',   'char', {'no', 'yes'});
ft_checkopt(cfg.mvpa, 'cv',         'char', {'kfold', 'leaveout', 'holdout','none'});

% flip dimensions such that trials come first
dat = dat';

% if cfg.dim has more than one non-singleton entry then the data has been
% multi-dimensional before. We must reshape it back to run it in MVPA-Light
n_extra_dim = sum(cfg.dim > 1);
if numel(cfg.dim) > 1
  dat = reshape(dat, [size(dat,1), cfg.dim]);
end

% design is the class labels (clabel)
clabel = design;

%% check for deprecated notation
if isfield(cfg, 'timextime')
    ft_warning('The parameter cfg.timextime has been removed, setting cfg.generalize = ''time'' instead')
    cfg.generalize = 'time';
end
if isfield(cfg, 'searchlight')
    ft_warning('The parameter cfg.searchlight has been removed, setting cfg.search = ''chan'' instead')
    cfg.search = {'chan'};
end

%% check if any cfg parameters are incompatible
if ~strcmp(cfg.generalize,'no') && isempty(strfind(cfg.dimord, cfg.generalize))
    ft_error('cfg.generalize has been set to ''%s'' but there is no such dimension since dimord is ''%s''', cfg.generalize, cfg.dimord)
end

if ~strcmp(cfg.generalize,'no') && n_extra_dim < 1
  ft_warning('cfg.generalize = ''%s'' but there does not seem to be such a dimension in the data', cfg.generalize);
  cfg.timextime = 'no';
end

%% extract dimensions from cfg.dimord or guess them
if isfield(cfg,'dimord')
    dimord_cell = [{'rpt'} strsplit(cfg.dimord,'_')];
else
    ft_warning('no dimord found, trying to guess it')
    % identify the calling function to know whether its timelock or freq data
    [calling, I] = dbstack;
    calling = calling(I+1).name;
    is_timelock = strcmp(calling, 'ft_timelockstatistics');
    is_freq     = strcmp(calling, 'ft_freqstatistics');
    dimord_cell = [];
    if is_timelock, dimord_cell= {'rpt' 'chan' 'time'};
    elseif is_freq, dimord_cell = {'rpt' 'chan' 'freq' 'time'}; end
    dimord_cell = dimord_cell(1:n_extra_dim+1);
end

% if data has been averaged over time/chans we must account for this in dimord
if isfield(cfg,'avgoverchan') && strcmp(cfg.avgoverchan,'yes'), dimord_cell(ismember(dimord_cell,'chan')) = []; end
if isfield(cfg,'avgovertime') && strcmp(cfg.avgovertime,'yes'), dimord_cell(ismember(dimord_cell,'time')) = []; end
if isfield(cfg,'avgoverfreq') && strcmp(cfg.avgoverfreq,'yes'), dimord_cell(ismember(dimord_cell,'freq')) = []; end

cfg.mvpa.dimension_names = get_nice_dimnames(dimord_cell); % nice names for console output in MVPA-Light
if numel(cfg.mvpa.dimension_names)==1
    cfg.mvpa.dimension_names = [cfg.mvpa.dimension_names {''}]; % fix to get correct feedback
end

%% define search dimensions
ft_checkopt(cfg,      'search',     {'char' 'cell' 'empty'});
if any(ismember(dimord_cell,'time')), 	 cfg.search = ft_getopt(cfg, 'search','time', true);
elseif any(ismember(dimord_cell,'freq')), cfg.search = ft_getopt(cfg, 'search','freq', true);
else,                                cfg.search = ft_getopt(cfg, 'search','',     true);
end

if isempty(cfg.search) 
    n_search = 0;
elseif ischar(cfg.search)
    n_search = 1;
    cfg.search = {cfg.search};
else
    n_search = numel(cfg.search);
end

if ~all(ismember(cfg.search, {'time' 'freq' 'chan' ''}))
    ft_error('cfg.search can only take the values: time freq chan')
end

dimord = strjoin(dimord_cell,'_');

%% neighbours
has_neighbours = ~isempty(cfg.mvpa.neighbours);

if has_neighbours
    % convert neighbours from struct array to matrix
    if isstruct(cfg.mvpa.neighbours)
        cfg.mvpa.neighbours = channelconnectivity(struct('neighbours',cfg.mvpa.neighbours, 'channel', {cfg.channel}));
    elseif iscell(cfg.mvpa.neighbours)
        for ii=1:numel(cfg.mvpa.neighbours)
            if isstruct(cfg.mvpa.neighbours{ii})
                cfg.mvpa.neighbours{ii} = channelconnectivity(struct('neighbours',cfg.mvpa.neighbours{ii}, 'channel', {cfg.channel}));
            end
        end
    end
end

%% Call MVPA-Light
% MVPA-Light has a number of different high-level functions. Here, we need
% to map the FieldTrip data to the most suited high-level function by
% looking at how many data dimensions there are, how many search dimensions
% are specified, and whether generalization is required.
label = [];
dim   = [];
stat_dimord = '';

if n_search == 0   
    % --- there are no search dimensions, so perform just a single cross-validation ---
    if n_extra_dim <= 1
        [perf, result] = mv_crossvalidate(cfg.mvpa, dat, clabel);
    else
        % we have to use mv_classify because there's multiple feature
        % dimensions that need to be flattened
        cfg.mvpa.sample_dimension  = 1;
        cfg.mvpa.feature_dimension = 2:ndims(dat); % all data dimensions are used as features
        [perf, result] = mv_classify(cfg.mvpa, dat, clabel);
    end

elseif strcmp(cfg.search{1}, 'chan') && (strcmp(dimord,'rpt_chan') || strcmp(dimord,'rpt_chan_time') || strcmp(dimord,'rpt_chan_time'))
        % --- searchlight across channels ---
        
        [perf, result] = mv_searchlight(cfg.mvpa, dat, clabel);
        
        % this preserves any spatial dimension, so no adjustment is done to a
        % channel list, if present
        if isfield(cfg, 'channel')
            label = cfg.channel;
        end
        
        if isfield(cfg, 'dim')
            dim = cfg.dim;
        end
        
        stat_dimord = 'chan';

elseif n_search==1 && ~has_neighbours && (strcmp(dimord,'rpt_chan_freq') || strcmp(dimord,'rpt_chan_time'))
    % --- 1 search dimension ---
    % Perform classification (or generalization) across time (or freq)
    
    if strcmp(cfg.generalize, 'time')
        % --- time x time generalization ---
        [perf, result] = mv_classify_timextime(cfg.mvpa, dat, clabel);
        
        % this does not preserve any spatial dimension, so label should be
        % adjusted
        label = squeezelabel(label, cfg);
        dim   = squeezedim(dim, cfg);
        stat_dimord = 'time_time';
        
    elseif strcmp(cfg.generalize, 'freq')
        
        error('TODO')
        stat_dimord = 'freq_freq';
        
    else
        % classification across time or frequency, but we can use
        % mv_classify_across_time in any case
        [perf, result] = mv_classify_across_time(cfg.mvpa, dat, clabel);
        
        % this does not preserve any spatial dimension, so label should be
        % adjusted
        label = squeezelabel(label, cfg);
        dim   = squeezedim(dim, cfg);
        stat_dimord = cfg.search{1};
        
%     else
%         % --- data has no time dimension, perform only one cross-validation ---
%         [perf, result] = mv_crossvalidate(cfg.mvpa, dat, clabel);
%         
%         % this does not preserve any spatial dimension, so label should be
%         % adjusted
%         label = squeezelabel(label, cfg);
%         dim   = squeezedim(dim, cfg);
%         stat_dimord = '';
    end
    
else
    % For all other cases, we use the general-purpose classification
    % function mv_classify. Just need to set up the dimensions correctly
    % according to dimord.
    %% TODO
    cfg.mvpa.sample_dimensions = 1;
    
    if ~strcmp(cfg.generalize,'no')
    end
    
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

if ~isempty(stat_dimord)
  stat.dimord = stat_dimord;
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

function dim_names_out = get_nice_dimnames(dim_names_in)
% Creates nicer dimnames 

dim_names_out = cell(numel(dim_names_in),1);
for ix=1:numel(dim_names_in)
    switch(dim_names_in{ix})
        case {'rpt','trial'}, dim_names_out{ix} = 'trials';
        case 'chan', dim_names_out{ix} = 'channels';
        case 'time', dim_names_out{ix} = 'time points';
        case 'freq', dim_names_out{ix} = 'frequencies';
        otherwise, dim_names_out{ix} = dim_names_in{ix};
    end
end

