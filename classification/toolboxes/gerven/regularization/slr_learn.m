function [model,diagnostics] = slr_learn(cfg,data)
% SLR_LEARN learns a sparse multinomial logistic regression model
% using L1/LP regularization
%
% Usage:
%
%   [model,diagnostics] = slr_learn(cfg,data)
%
%   learns a logistic regression model that is regularized using L1/LP
%   regularization. 
%
%   model is the weight matrix of the best model in case cfg.folds <= 1
%
% Notes: 
%   - data should be centered for good performance
%   - cfg.lambda = 0 tries to perform unregularized logistic regression 
%       using the statistics toolbox 
%
% SEE ALSO:
%   regularize.m
%
% Copyright (c) 2008, Marcel van Gerven
% F.C. Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% $Log: slr_learn.m,v $
% Revision 1.8  2008/05/19 09:36:45  marvger
% regular update
%
% Revision 1.7  2008/04/16 09:17:49  marvger
% large update
%
% Revision 1.6  2008/04/09 07:36:09  marvger
% regular update
%
% Revision 1.5  2008/03/17 12:47:41  marvger
% regular update
%
% Revision 1.4  2008/03/03 16:32:21  marvger
% changed lr_reg_grad computation exception: partial reason for problems
% when p becomes large
%
% Revision 1.3  2008/03/01 16:38:28  marvger
% removed somes files
%
% Revision 1.2  2008/02/29 12:28:39  marvger
% added unregularized logistic regression
%
% Revision 1.1.1.1  2008/02/27 14:42:55  roboos
% Marcel van Gerven, version 27 Feb 2008
%
%

if nargin < 2, error('cfg and data must be specified'); end

% for multitask learning, we call a specialized algorithm
if iscell(data)
   error('dont know how to handle cell data; call slr_learn_transfer instead');
end

if ~isfield(cfg,'cvar'), cfg.cvar = 1; end

if ~isfield(cfg,'metric'), cfg.metric = 'accuracy'; end

if ~isfield(cfg,'folds'), cfg.folds = 0.8; end

if cfg.folds == 1

    testdata = data;
    traindata = data;
elseif ~isscalar(cfg.folds) % explicit test set
    
    testdata = cfg.folds;
    traindata = data;

else
    % choose random permutation for training and testing
    cfg.nexamples = size(data,1);
    perm = randperm(cfg.nexamples);
    testdata = data(perm((floor(cfg.folds * cfg.nexamples)+1):cfg.nexamples),:);
    traindata = data(perm(1:floor(cfg.folds * cfg.nexamples)),:);
end

% learn
if isfield(cfg,'lambda') && length(cfg.lambda) == 1 && cfg.lambda == 0

    % proprietary code
    if exist('minimize','file')
                
        cfg.nexamples = size(traindata,1);
        targets = traindata(:,cfg.cvar);
        traindata = [traindata(:,[1:(cfg.cvar-1) (cfg.cvar+1):end]) ones(cfg.nexamples,1)];
        cfg.nclasses = max(targets);
        cfg.nfeatures = size(traindata,2);
        classidxs = (1:cfg.nexamples)' + (targets - 1) .* cfg.nexamples;
        ptargets = zeros(cfg.nexamples,cfg.nclasses);
        ptargets(classidxs) = 1;
        targets = (1:cfg.nexamples)' + (targets - 1) * cfg.nexamples;

        w = zeros(cfg.nclasses,cfg.nfeatures); w = w(:);
        [w, fw, i] = minimize(w,'logreg',50,traindata,targets,ptargets,cfg.nclasses);

        diagnostics.w = w;
        diagnostics.fw = fw;
        diagnostics.i = i;
        
        res = { reshape(w,cfg.nclasses,cfg.nfeatures); };
        
    else

        % try unregularized logistic regression
        if license('test','statistics_toolbox') && size(data,2) <= 300

            % try matlab native code
            [B,dev,diagnostics] = mnrfit(traindata(:,[1:(cfg.cvar-1) (cfg.cvar+1):end]),traindata(:,cfg.cvar));
            diagnostics.dev = dev;

            % construct model
            res = { [transpose([B(2:end,:); B(1,:)]); zeros(1,size(B,1))] };

        else % slow steepest descent

            [res,diagnostics] = regularize_lr(cfg,traindata);
        end
    end
else
    [res,diagnostics] = regularize_lr(cfg,traindata);    
end

% keep accuracy results
theaccuracy = zeros(1,length(res));
trainaccuracy = zeros(1,length(res));

% compute accuracies
for r=1:length(res)

    % order from small lambda to large lambda
    res2 = res{r};

    % test
    p2 = slr_classify([testdata(:,[1:(cfg.cvar-1) (cfg.cvar+1):end]) ones(size(testdata,1),1)],res2);

    % test on training data
    p3 = slr_classify([traindata(:,[1:(cfg.cvar-1) (cfg.cvar+1):end]) ones(size(traindata,1),1)],res2);
    
    if strcmp(cfg.metric,'loglik')

        % use log loss as the objective function
        td = testdata(:,cfg.cvar);
        for k=1:length(td)
            theaccuracy(r) = theaccuracy(r) + log(p2(k,td(k)));
        end
        td = traindata(:,cfg.cvar);
        for k=1:length(td)
            trainaccuracy(r) = trainaccuracy(r) + log(p3(k,td(k)));
        end
    else
        % use accuracy as objective function
        [m,i] = max(p2,[],2);
        theaccuracy(r) = sum((i == testdata(:,cfg.cvar))) / length(testdata(:,cfg.cvar));
        [m,i] = max(p3,[],2);
        trainaccuracy(r) = sum((i == traindata(:,cfg.cvar))) / length(traindata(:,cfg.cvar));
    end
    
end

% keep training performance for evaluation purpose
diagnostics.trainperformance = trainaccuracy;

% in this case it is useful to retain raw results for
% further processing
diagnostics.path = res;
diagnostics.testperformance = theaccuracy;

[m,idx] = max(theaccuracy);

if ~isempty(idx)
    model = diagnostics.path{idx}; % retrieve the weight matrix of the best result
else
    error('no suitable model found');
end

