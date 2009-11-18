function [model,diagnostics] = slr_learn_transfer(cfg,data)
% SLR_LEARN_TRANSFER transfer learns a sparse multinomial logistic regression model
% using L1/LP regularization
%
% Usage:
%
%   [model,diagnostics] = slr_learn_transfer(cfg,data)
%
%   transfer learns a logistic regression model that is regularized using L1/LP
%   regularization. 
%
%   model is the weight matrix of the best model in case cfg.folds <= 1
%
%   cfg.soft = false chooses between standard or soft transfer learning
%
% Note: data should be centered for good performance
% 
% Copyright (c) 2008, Marcel van Gerven
% F.C. Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% $Log: slr_learn_transfer.m,v $
% Revision 1.2  2008/03/01 16:38:28  marvger
% removed somes files
%
% Revision 1.1.1.1  2008/02/27 14:42:55  roboos
% Marcel van Gerven, version 27 Feb 2008
%
%

if ~isfield(cfg,'soft'), cfg.soft = false; end

if nargin < 2, error('cfg and data must be specified'); end

if ~iscell(data), error('expecting data as a cell array for each task'); end

if ~isfield(cfg,'cvar'), cfg.cvar = 1; end

if ~isfield(cfg,'metric'), cfg.metric = 'accuracy'; end

if ~isfield(cfg,'folds'), cfg.folds = 0.8; end

cfg.ntasks = length(data);
cfg.nexamples = cell(1,cfg.ntasks);
traindata = cell(1,cfg.ntasks);
testdata = cell(1,cfg.ntasks);

% make random permutations for train and test data
for s=1:cfg.ntasks    
    if isscalar(cfg.folds)
        cfg.nexamples{s} = size(data{s},1);
        perm = randperm(cfg.nexamples{s});
        testdata{s} = data{s}(perm((floor(cfg.folds * cfg.nexamples{s})+1):cfg.nexamples{s}),:);
        traindata{s} = data{s}(perm(1:floor(cfg.folds * cfg.nexamples{s})),:);
    else
        % explicit test set
        testdata = cfg.folds;
        traindata = data;
    end
end

% learn
if cfg.soft
    [res,diagnostics] = regularize_lr_soft_transfer(cfg,traindata);
else
    [res,diagnostics] = regularize_lr_transfer(cfg,traindata);
end

% keep accuracy results
theaccuracy = zeros(1,length(res));
trainaccuracy = zeros(1,length(res));

% compute accuracies
for r=1:length(res)

    % order from small lambda to large lambda
    res2 = res{r};

    % the accuracy is averaged over the tasks
    theaccuracy(r) = 0;
    trainaccuracy(r) = 0;
    for s=1:cfg.ntasks
        
        % test
        p2 = slr_classify([testdata{s}(:,[1:(cfg.cvar-1) (cfg.cvar+1):end]) ones(size(testdata{s},1),1)],squeeze(res2(s,:,:)));

        % test on training data
        p3 = slr_classify([traindata{s}(:,[1:(cfg.cvar-1) (cfg.cvar+1):end]) ones(size(traindata{s},1),1)],squeeze(res2(s,:,:)));
        
        if strcmp(cfg.metric,'loglik')
            % use log loss as the objective function
            td = testdata{s}(:,cfg.cvar);
            for k=1:length(td)
                theaccuracy(r) = theaccuracy(r) + log(p2(k,td(k)));
            end
            td = traindata{s}(:,cfg.cvar);
            for k=1:length(td)
                trainaccuracy(r) = trainaccuracy(r) + log(p3(k,td(k)));
            end
        else
            % use accuracy as objective function
            [m,i] = max(p2,[],2);
            theaccuracy(r) = theaccuracy(r) + sum((i == testdata{s}(:,cfg.cvar))) / length(testdata{s}(:,cfg.cvar));
            [m,i] = max(p3,[],2);
            trainaccuracy(r) = trainaccuracy(r) + sum((i == traindata{s}(:,cfg.cvar))) / length(traindata{s}(:,cfg.cvar));
        end
    
    end
    theaccuracy(r) = theaccuracy(r)/cfg.ntasks;
    trainaccuracy(r) = trainaccuracy(r)/cfg.ntasks;
end

% keep training performance for evaluation purpose
diagnostics.trainperformance = trainaccuracy;

% in this case it is useful to retain raw results for
% further processing
diagnostics.path = res;
diagnostics.testperformance = theaccuracy;

[m,idx] = max(theaccuracy);

if ~isempty(idx)
    
    % retrieve the weight matrix of the best result
    tmpmodel = diagnostics.path{idx}; 
    
    % transform to a result per subject
    model = cell(1,cfg.ntasks);
    for s=1:cfg.ntasks
        model{s} = squeeze(tmpmodel(s,:,:));
    end
    
else
    error('no suitable model found');
end
