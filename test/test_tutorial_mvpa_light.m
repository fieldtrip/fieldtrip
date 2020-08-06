function test_tutorial_mvpa_light

% MEM 9gb
% WALLTIME 00:20:00
% DEPENDENCY ft_timelockstatisitcs ft_statistics_mvpa

%
%% Classification of event related MEG data using MVPA-Light
%
%% # Introduction
%
% The objective of this tutorial is to give an introduction to the classification of event related
% data using the [MVPA-Light](https://github.com/treder/MVPA-Light) toolbox. For a general introduction and background on multivariate analysis, refer to the
% and the [MVPA-Light readme file](https://github.com/treder/MVPA-Light/blob/master/README.md).
% This tutorial builds on skills acquired in the [preprocessing](/tutorial/preprocessing), [event related averaging](/tutorial/eventrelatedaveraging) and [time-frequency analysis](/tutorial/timefrequencyanalysis) tutorials.
%
%
%% # Installation
%
% MVPA-Light is a stand-alone Matlab toolbox for multivariate pattern analysis. FieldTrip provides a high-level interface to its functions so one does not need to directly interact with the toolbox.
% However, it needs to be installed and included in
% Matlab's search path. To this end, [follow the installation instructions](https://github.com/treder/MVPA-Light#installation-) on its Github
% page.
%
%% # Procedure
%
% We will use classifiers to analyze the [MEG-language dataset](/faq/what_types_of_datasets_and_their_respective_analyses_are_used_on_fieldtrip) which
% features one subject with three types of trials: fully incongruent (FIC), fully congruent (FC), and
% initially congruent (IC). These three classes are stored in different files available here:
%
% The data can be loaded into MATLAB using
%
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataFIC_LP'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataFC_LP'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataIC_LP'));

% To get started, we investigate whether we can discriminate between the three classes
% FIC, FC, and IC, using the average activity in the 0.5-0.7 s post-stimulus interval.
%
% After this, we will focus on two out of the three classes, namely FIC vs FC, and we will investigate the following questions:
%
%* At *what times* ('when') in a trial can one discriminate between FIC and FC?
%* At *which sensor locations* ('where') can one discriminate between FIC and FC?
%* Which representations that discriminate between FIC and FC [*generalise across time*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5635958/)?
%
% Note that the classification is performed for a single subject using single trials.
%
%% # Classification in the 0.5-0.7 s interval
%
% We will use `ft_timelockstatistics` to determine the classification accuracy between the three classes FIC, FC, and IC. As features, the average activity in each MEG channel in the the 0.5-0.7 s interval is used. In each trial, this yields 149 features for the classifier, one feature per MEG channel. Let us first determine the number of trials in each class:
%
nFIC = numel(dataFIC_LP.trial);
nFC = numel(dataFC_LP.trial);
nIC = numel(dataIC_LP.trial);

% Define the configuration struct
%
cfg = [] ;
cfg.method      = 'mvpa';
cfg.mvpa.classifier  = 'multiclass_lda';
cfg.mvpa.metric      = 'accuracy';
cfg.mvpa.k           = 3;
cfg.latency     = [0.5, 0.7];
cfg.avgovertime = 'yes';
cfg.design      = [ones(nFIC,1); 2*ones(nFC,1); 3*ones(nIC,1)];

% Let us unpack this:
%
%* `cfg.mvpa.classifier` indicates which classifier we want to use. Here, we use multi-class Linear Discriminant Analysis (LDA).  [Click here](https://github.com/treder/MVPA-Light#classifiers-for-two-classes) for a full list of available classifiers.
%* `cfg.metric` indicates the metric we use to measure classifier performance. Here, *classification accuracy* is used. Other metrics such as AUC and F1-score are available. [Click here](https://github.com/treder/MVPA-Light#classifier-performance-metrics) for a full list of available metrics.
%* `cfg.mvpa.k` specifies the number of folds used to calculate the cross-validated performance. Cross-validation is explained in more detail in the next section.
%* `cfg.latency` restricts the classification analysis to a specific time window (here 0.5-0.7s).
%* `cfg.avgovertime` specifies whether the activity in latency window should be averaged prior to classification. If `'no'`, a separate classification is performed for every time point (see section *Classification across time*).
%* `cfg.design` specifies the vector of *class labels*. Class labels indicate which class (or experimental condition) trials belong to. The task of the classifier is to predict these class labels given the data. To this end, we create a vector with *1*'s
% for the trials belonging to class 1, *2*'s for trials
% belonging to class 2, and so on. For the [MEG-language dataset](/faq/what_types_of_datasets_and_their_respective_analyses_are_used_on_fieldtrip),
% there is three classes, namely FIC (class 1), FC (class 2), and IC (class 3).
%
% Now call
%
stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP, dataIC_LP);

% to perform the classification analysis. It is important to make sure that the order of class labels (FIC, FC, IC) matches the order that the datasets are passed in to `ft_timelockstatistics`. Let us print the resulting classification
% accuracy
%
fprintf('Classification accuracy: %0.2f\n', stat.accuracy)

% For multi-class problems, the [confusion matrix](https://en.wikipedia.org/wiki/Confusion_matrix) is
% a useful metric. In a confusion matrix, rows correspond to the true class labels,
% columns correspond to predicted class labels. The (i,j)-th element gives the
% proportion of samples of class i that have been classified as class j. Consequently,
% the diagonal of the confusion matrix contains the proportion of correct classifications. Off-diagonal elements specify the misclassifications. To obtain
% the confusion matrix, all we need to do is to change the |metric| field:
%
cfg.mvpa.metric      = 'confusion';
stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP, dataIC_LP);

stat.confusion

% Looking at the diagonal of the matrix tells us that the classifier is doing better
% at predicting classes 1 and 2 than it is at correctly predicting class 3.
% For a simple visualisation of this result, we can use a plotting function in  [MVPA-Light](https://github.com/treder/MVPA-Light)
% called [`mv_plot_result`](https://github.com/treder/MVPA-Light/blob/master/plot/mv_plot_result.m).
% It takes the result structure returned in `stat.mvpa_result` which contains the
% classification results in a format required by the function.
%
mv_plot_result(stat.mvpa)

%
%
%% ## Cross-validation
%
% To obtain a realistic estimate of classifier performance and control for overfitting, a classifier should be tested on an independent dataset that has not been used for training. In most neuroimaging experiments, there is only one dataset with a restricted number of trials. K-fold [cross-validation](https://en.wikipedia.org/wiki/Cross-validation) makes efficient use of this data by splitting it into k different folds. In each iteration, one of the k folds is held out and used as test set, whereas all other folds are used for training the model. This process is repeated until every fold has been used as test set once. Cross-validation is controlled by the following parameters:
%
%* `cfg.mvpa.cv`: cross-validation type, either 'kfold', 'leaveout' or 'holdout' (default 'kfold')
%* `cfg.mvpa.k`: number of folds or partitions in k-fold cross-validation (default 5)
%* `cfg.mvpa.repeat`: number of times the whole cross-validation analysis is repeated with new randomly assigned folds (default 5)
%* `cfg.mvpa.p`: if `cfg.mvpa.cv` is 'holdout', |p| is the fraction of test samples (default 0.1)
%* `cfg.mvpa.stratify`: if 1, the class proportions are approximately preserved in each test fold (default 1)
%
% The total number of training and testing iterations is equal to `cfg.k * cfg.repeat`. The result returned by `ft_timelockstatistics` is an average
% across the test folds.
%
%% ### Exercise 1
%
% What is the effect of setting k to a very large vs very small value? Why is it
% useful to repeat the cross-validation multiple times?
%
%% # Classification across time ('when')
%
% Many neuroimaging datasets have a 3-D structure (trials x channels x time). Classification across time can help identify the time points in a trial *when* discriminative information shows up. To this end, classification is performed for each time point separately. First, we need to make sure that the time dimension is not averaged out. We can set `cfg.avgovertime = 'no'`, but since the default value is `'no'` we can simply omit this parameter.
%
cfg = [] ;  
cfg.method           = 'mvpa';
cfg.mvpa.classifier  = 'lda';
cfg.mvpa.metric      = 'auc';
cfg.mvpa.k           = 10;
cfg.mvpa.repeat      = 2;
cfg.design           = [ones(nFIC,1); 2*ones(nFC,1)];

% For simplicity, we will limit ourselves to comparing only FIC and FC. As classifier,
% we use Linear Discriminant Analysis (LDA). As metric, we use area under the ROC curve (AUC).
% It is calculated using 10-fold cross-validation with 2 repetitions.
%
stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP);

% Note that the metric is now a vector with 900 values, one AUC value for each time point in the trial.
% It can be plotted as a function of time using
%
figure; plot(stat.auc)

% For a slightly nicer plot, one can again use `mv_plot_result`. As an additional parameter,
% we can pass the values for the time axis. This makes sure that the x-axis is formatted
% correctly.
%
mv_plot_result(stat.mvpa, dataFC_LP.time{1})

% The resultant plot shows AUC across time in the trial. The shaded area
% is the standard deviation of the AUC metric across the different test sets in the
% cross-validation.
%
%
%
%
%
%% ### Exercise 2
%
% Perform classification across time using all three classes FIC, FC, and IC. As
% classifier, use kernel FDA. As metric, use classification accuracy.
%
%
%% # Searchlight analysis ('where')
%
% Which channels contribute most to classification performance? The answer to this question can be used to better interpret the data or to perform feature selection. To this end, we will perform classification for each feature separately. The result of the searchlight analysis is a classification performance measure for each channel.
%
cfg = [] ;  
cfg.method      = 'mvpa';
cfg.searchlight = 'yes';
cfg.design      = [ones(nFIC,1); 2*ones(nFC,1)];
cfg.latency     = [0.3, 0.7];
cfg.avgovertime = 'yes';
stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP);

% Since we did not specify a classifier and a metric, the default values (LDA as classifier and classification accuracy as metric)
% are used. In searchlight analysis, the *time points* in a trial are used as
% features, for each channel separately. Set `cfg.latency` to restrict the analysis to
% a specific time window. Set `cfg.avgovertime='yes'` if you prefer to average the values in the time window to a single feature.
%
% Since a classification result is obtained for each channel, classification accuracy can be plotted as a topography. 
% Then call `ft_topoplotER` to do the plotting.
%

cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = 'CTF151_helmet.mat';            
cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
ft_topoplotER(cfg, stat);

%
%
% In the previous analysis, classification has been performed for each channel individually.
% However, since the MEG channels have a spatial structure,
% one can also consider groups of neighbouring channels in the searchlight. To do this, we must provide
% a distance matrix that specifies which channels are neighbours of each other.
%
%%% Get layout
cfg = [];
cfg.layout      = 'CTF151_helmet.mat';
cfg.skipscale   = 'yes';
cfg.skipcomnt   = 'yes';
cfg.channel     = dataFIC_LP.label;
lay = ft_prepare_layout(cfg);

%%% Get distance matrix
nb_mat = squareform(pdist(lay.pos));

% We are now ready to re-run the searchlight analysis. We pass the neighbourhood distance matrix
% via the parameter `cfg.nb`. By setting `cfg.size = 3`  in every iteration
% the target channel is considered together with its 3 closest neighbouring channels.
%
  cfg = [] ;  
  cfg.method      = 'mvpa';
  cfg.searchlight = 'yes';
  cfg.design      = [ones(nFIC,1); 2*ones(nFC,1)];
  cfg.latency     = [0.3, 0.7];
  cfg.avgovertime = 'yes';

  cfg.mvpa.nb          = nb_mat;
  cfg.mvpa.size        = 3;

  stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP);

  cfg              = [];
  cfg.parameter    = 'accuracy';
  cfg.layout       = 'CTF151_helmet.mat';            
  cfg.xlim         = [0, 0];
  cfg.colorbar     = 'yes';
  ft_topoplotER(cfg, stat);

  
% As expected, the resultant topography is slightly more smeared out. Peak classification accuracy is higher which is due to the classifier now combining information across neighbouring channels.
%
%
%
%% # Time generalisation (time x time classification)
%
% Classification across time does not give insight into whether information is shared across different time points. For example, is the information that the classifier uses early in a trial (t=80 ms) the same that it uses later (t=300ms)? In time generalisation, this question is answered by training the classifier at a certain time point t. The classifer is then tested at the same time point t but it is also tested at all other time points in the trial ([King and Dehaene, 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5635958/)). To perform
% time x time classification, we only need to set the `cfg.timextime` parameter: 
%
%
cfg = [] ;  
cfg.method      = 'mvpa';
cfg.timextime   = 'yes';
cfg.design      = [ones(nFIC,1); 2*ones(nFC,1)];

stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP);

% It returns a 2-D matrix of classification performance, with performance calculated for each combination of training time point and testing time point. We plot the
% result using [`mv_plot_result`](https://github.com/treder/MVPA-Light/blob/master/plot/mv_plot_result.m). As parameters, we pass the classification result and two additional parameters specifying the x-axis and y-axis.
%
mv_plot_result(stat.mvpa, dataFC_LP.time{1}, dataFC_LP.time{1})

% In the resultant plot, each row (corresponding to a value of on the y-axis)  corresponds to the
% time point at which the classifier was trained. Each point on the x-axis corresponds
% to a time point at which the respective classifier was tested. Clearly, the classifier attains peak performance roughly in the 0.45-0.65s period.
%
%
%
%
%
%% # Advanced topics
%
% In this section we address slightly more advanced topics that might become important
% once one starts using MVPA on a regular basis.
%
%% ## Hyperparameters
%
% Many classifiers have parameters that control their properties and need to
% be set by the user, so-called *hyperparameters*. For a list of hyperparameters for each classifier, see the respective train_
% functions in the [classifier folder](https://github.com/treder/MVPA-Light/tree/master/classifier).
% Hyperparameters can be set using the |param| substruct. For instance, in Support
% Vector Machines (SVM) the kernel is a hyperparameter and |gamma| controls the
% kernel width for an RBF kernel.
%
%
cfg.mvpa.param           = [];
cfg.mvpa.param.kernel    = 'rbf';
cfg.mvpa.param.gamma     = 1;

% See [train_svm](https://github.com/treder/MVPA-Light/blob/master/classifier/train_svm.m)) for a list of SVM hyperparameters and their default values.
% To give another example, in LDA the |lambda| parameter controls the amount of regularisation of the covariance matrix.
%
cfg.param           = [];
cfg.param.lambda    = 'auto';

% See [train_lda](https://github.com/treder/MVPA-Light/blob/master/classifier/train_svm.m)) for a list of LDA hyperparameters and their default values.
% In many cases the default values suffice.
%
%
%% ## Unbalanced classes
%
% Classes are unbalanced when one class contains more instances than another class.
% Unbalanced classes can distort some of the classification metrics. For instance,
% if 90\% of the metrics
%
%
%% ## Classifier weights vs activation patterns
%
% TODO
%
%
%
%% # Summary
%
% In this tutorial, classification across time, searchlight analysis, and time generalisation
% information within MEG trials. It is worth stressing that searchlight analysis and classification across time yield complementary information: in searchlight analysis,
% the time points serve as features and classification is performed for each channel separately.
% In classification across time, the channels serve as features and classification is performed for each time point separately.
%
% TODO:
% advanced topics
% statistical analysis
%
