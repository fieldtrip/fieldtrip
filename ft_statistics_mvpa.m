function [stat, cfg] = ft_statistics_mvpa(cfg, dat, design)

% FT_STATISTICS_MVPA performs multivariate pattern classification 
% on the data. If the data has not been averaged over time, classification
% is performed separately for every time point. Additionally, searchlight
% analysis can be performed (classification for each channel/voxel
% separately), or time x time generalisation.
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
%   cfg.classifier        = string, classifier to use
%                    'lda'          Regularised linear discriminant analysis
%                                   (LDA) (for two classes)
%                    'multiclass_lda' LDA for more than two classes
%                    'logreg'       Logistic regression with L2-regularisation
%                    'svm'          Support Vector Machine (SVM) with L2-regularisation
%                    'ensemble'     Ensemble of classifiers. Any of the other
%                                   classifiers can be used as a learner.
%   cfg.param             = struct, structure with hyperparameters for the 
%                           classifier
%   cfg.statistic         = string, performance statistic (called 'metric' 
%                           within MVPA-Light). Possible statistics are:
%                           'accuracy'
%                           'auc'
%
%   cfg.searchlight       = 'yes' or 'no', performs searchlight analysis
%                           (default 'no')
%   cfg.timextime         = 'yes' or 'no', performs time x time
%                           generalisation. In other words, the classifier
%                           is trained at each time point and tested at
%                           every time point. The result is a time x time
%                           matrix of classification performance.
%                           (default 'no')
%                           Note that searchlight and timextime cannot be
%                           run simultaneously (at least one option needs
%                           to be set to 'no').
%
% .balance      - for imbalanced data with a minority and a majority class.
%                 'oversample' oversamples the minority class
%                 'undersample' undersamples the minority class
%                 such that both classes have the same number of samples
%                 (default 'none'). Note that for we undersample at the
%                 level of the repeats, whereas we oversample within each
%                 training set (for an explanation see mv_balance_classes).
%                 You can also give an integer number for undersampling.
%                 The samples will be reduced to this number. Note that
%                 concurrent over/undersampling (oversampling of the
%                 smaller class, undersampling of the larger class) is not
%                 supported at the moment
% .replace      - if balance is set to 'oversample' or 'undersample',
%                 replace deteremines whether data is drawn with
%                 replacement (default 1)
% .normalise    - normalises the data across samples, for each time point 
%                 and each feature separately, using 'zscore' or 'demean' 
%                 (default 'zscore'). Set to 'none' or [] to avoid normalisation.
% .feedback     - print feedback on the console (default 1)
%
% To obtain a realistic estimate of classification performance,
% cross-validation is used. It is controlled by the following parameters:
%   cfg.cv              cross-validation type, either 'kfold', 'leaveout' 
%                       or 'holdout'. If 'none', no cross-validation is
%                       used and the classifier is tested on the training
%                       set. (default 'kfold')
%   cfg.k               number of folds in k-fold cross-validation (default 5)
%   cfg.repeat          number of times the cross-validation is repeated 
%                       with new randomly assigned folds (default 5)
%   cfg.p               if cfg.cv is 'holdout', p is the fraction of test 
%                       samples (default 0.1)
%   cfg.stratify        if 1, the class proportions are approximately 
%                       preserved in each test fold (default 1)
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
% additional SEARCHLIGHT ANALYSIS parameters:
% .nb          - [features x features] matrix specifying which features
%                are neighbours of each other.
%                          - EITHER - 
%                a GRAPH consisting of 0's and 1's. A 1 in the 
%                (i,j)-th element signifies that feature i and feature j 
%                are neighbours, and a 0 means they are not neighbours
%                            - OR -
%                a DISTANCE MATRIX, where larger values mean larger distance.
%                If no matrix is provided, every feature is only neighbour
%                to itself and classification is performed for each feature 
%                separately.
% .size        - if a nb matrix is provided, size defines the 
%                size of the 'neighbourhood' of a feature.
%                if nb is a graph, it gives the number of steps taken 
%                     through the nb matrix to find neighbours:
%                     0: only the feature itself is considered (no neighbours)
%                     1: the feature and its immediate neighbours
%                     2: the feature, its neighbours, and its neighbours'
%                     neighbours
%                     3+: neighbours of neighbours of neighbours etc
%                     (default 1)
%                if nb is a distance matrix, size defines the number of
%                     neighbouring features that enter the classification
%                     0: only the feature itself is considered (no neighbours)
%                     1: the feature and its first closest neighbour 
%                        according to the distance matrix
%                     2+: the 2 closest neighbours etc.
%
% TIME x TIME GENERALISATION:
% -- TODO: how to specifiy the second dataset --
%
% Returns:
%   stat.statistic    = the statistics to report


%   stat.model        = the models associated with this multivariate analysis

% do a sanity check on the input data
assert(isnumeric(dat),    'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');
assert(isnumeric(design), 'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');

% flip dimensions such that the number of trials comes first
dat = dat';

% if cfg.dim has two entries which are non-singleton then the data has been 
% 3D before. We must reshape it from 2D to 3D to run it in MVPA-Light
data_is_3D = (numel(cfg.dim) > 1 && all(cfg.dim > 1));  % checks whether data was 3D before being reshaped
if data_is_3D
    dat = reshape(dat, size(dat,1), cfg.dim(1), cfg.dim(2));
end

% check if the input cfg is valid for this function
% cfg = ft_checkconfig(cfg, 'renamed',  {'hyperparameter', 'hyperparam'});

y = cfg.design;

%% build up config struct for MVPA-Light
mvcfg = keepfields(cfg, {'balance','replace','normalise', ...
                         'cv','k','repeat','p','stratify',...
                         'classifier','param','size','nb'});
mvcfg.metric          = ft_getopt(cfg, 'statistic','acc');
mvcfg.feedback        = ft_getopt(cfg, 'feedback','yes');

% cfg: set defaults
cfg.searchlight     = ft_getopt(cfg, 'searchlight','no');
cfg.timextime       = ft_getopt(cfg, 'timextime','no');

% % cfg: set defaults
% mvcfg.classifier      = ft_getopt(cfg, 'classifier','lda');
% mvcfg.param           = ft_getopt(cfg, 'param', []);
% mvcfg.metric          = ft_getopt(cfg, 'statistic','acc');
% mvcfg.feedback        = ft_getopt(cfg, 'feedback','yes');
% mvcfg.searchlight     = ft_getopt(cfg, 'searchlight','no');
% mvcfg.timextime       = ft_getopt(cfg, 'timextime','no');
% 
% % set cross-validation defaults
% mvcfg.cv              = ft_getopt(cfg, 'cv','kfold');
% mvcfg.k               = ft_getopt(cfg, 'k', 5);
% mvcfg.repeat          = ft_getopt(cfg, 'repeat', 5);
% mvcfg.p               = ft_getopt(cfg, 'p', 0.1);
% mvcfg.stratify        = ft_getopt(cfg, 'stratify', 1);

% translate parameter value into MVPA-Light notation [0 or 1]
if ischar(mvcfg.feedback)
    mvcfg.feedback = strcmp(mvcfg.feedback,'yes');
end

%% perform sanity checks on parameters
if strcmp(cfg.timextime,'yes') && strcmp(cfg.searchligh,'yes')
    ft_error('you should not set timextime = ''yes'' and searchlight = ''yes'' simultaneously')
end

% timextime = 'yes' but data is not 3D we should change timextime to 'no'
if strcmp(cfg.timextime,'yes') && ~data_is_3D 
    ft_warning('timextime = ''yes'' but data has no time dimension, setting timextime = ''no''');
    cfg.timextime = 'no';
end

%% Call MVPA-Light 
if strcmp(cfg.searchlight, 'yes')
    % --- searchlight analysis ---
    perf = mv_searchlight(mvcfg, dat, y);
    
elseif strcmp(cfg.timextime, 'yes')
    % --- time x time generalisation ---
    perf = mv_classify_timextime(mvcfg, dat, y);
    
elseif data_is_3D
    % --- classification across time ---
    perf = mv_classify_across_time(mvcfg, dat, y);

else
    % --- data has no time dimension, perform only cross-validation ---
    perf = mv_crossvalidate(mvcfg, dat, y);
end

%% setup stat struct
stat = [];
if ~iscell(cfg.statistic), cfg.statistic = {cfg.statistic}; end
if ~iscell(perf), perf = {perf}; end
for mm=1:numel(perf) 
    stat.statistic.(cfg.statistic{mm}) = perf{mm};
end