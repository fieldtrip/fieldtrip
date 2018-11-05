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
%   cfg.param             = struct, structure with hyperparameters for the 
%                           classifier
%   cfg.metric            = string, performance metric, possible metrics
%                           are:
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
%                           run simultanesouly (at least one option needs
%                           to be set to 'no').
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
% The different classifiers that are implemented are
%   cfg.classifier = 'lda'          Regularised linear discriminant analysis
%                                   (LDA) (for two classes)
%                    'multiclass_lda' LDA for more than two classes
%                    'logreg'       Logistic regression with L2-regularisation
%                    'svm'          Support Vector Machine (SVM) with L2-regularisation
%                    'ensemble'     Ensemble of classifiers. Any of the other
%                                   classifiers can be used as a learner.
%
% More information about each classifier is found in the documentation of
% MVPA-Light (github.com/treder/MVPA-Light/).
%
% HYPERPARAMETER:
% Each classifier comes with its own set of hyperparameters, such as 
% regularisation parameters and the kernel.
% The specification of the hyperparameters is found in the training function
% for each classifier at github.com/treder/MVPA-Light/tree/master/classifier
% If not specified, default values are used for the hyperparameters.
%
% SEARCHLIGHT ANALYSIS:
% neighbourhood matrix ...
%
% TIME x TIME GENERALISATION:
% -- how to specifiy the second dataset ?? --
%
% Returns:
%   stat.statistic    = the statistics to report


%   stat.model        = the models associated with this multivariate analysis

% do a sanity check on the input data
assert(isnumeric(dat),    'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');
assert(isnumeric(design), 'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');

% check if the input cfg is valid for this function
% cfg = ft_checkconfig(cfg, 'renamed',  {'hyperparameter', 'hyperparam'});

% set the defaults
cfg.classifier      = ft_getopt(cfg, 'classifier','lda');
cfg.param           = ft_getopt(cfg, 'param', []);
cfg.metric          = ft_getopt(cfg, 'metric','acc');
cfg.feedback        = ft_getopt(cfg, 'feedback','yes');
cfg.searchlight     = ft_getopt(cfg, 'searchlight','no');
cfg.timextime       = ft_getopt(cfg, 'timextime','no');

% set the cross-validation defaults
cfg.cv              = ft_getopt(cfg, 'cv','kfold');
cfg.k               = ft_getopt(cfg, 'k', 5);
cfg.repeat          = ft_getopt(cfg, 'repeat', 5);
cfg.p               = ft_getopt(cfg, 'p', 0.1);
cfg.stratify        = ft_getopt(cfg, 'stratify', 1);

% translate parameter value into MVPA-Light syntax [binary instead of string]
if ischar(cfg.feedback)
    cfg.feedback = strcmp(cfg.feedback,'yes');
end

%% Call MVPA-Light 
if strcmp(cfg.searchlight, 'yes')
    % --- searchlight analysis ---
    
    % if avgovertime , we have only one time feature ...
    
elseif strcmp(cfg.timextime, 'yes')
    % --- time x time generalisation ---
    
else
    % --- classification across time ---
end