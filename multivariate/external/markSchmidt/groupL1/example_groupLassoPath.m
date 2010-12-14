%% Demo of computing group Lasso regularization path with iterative
%% soft-thresholding with Barzilai-Borwein step sizes

clear all
close all

%% Generate Synthetic Data

% Generate categorical features
nInstances = 500;
nStates = [3 3 2 3 3 5 4 5 5 6 10 3 3 4 5 2 2 6 8 9 2 7]; % Number of discrete states for each categorical feature
X = zeros(nInstances,length(nStates));
offset = 0;
for i = 1:nInstances
    for s = 1:length(nStates)
        prob_s = rand(nStates(s),1);
        prob_s = prob_s/sum(prob_s);
        X(i,s) = sampleDiscrete(prob_s);
    end
end

% Make indicator variable encoding of categorical features
X_ind = zeros(nInstances,sum(nStates));
for s = 1:length(nStates)
    for i = 1:nInstances
        X_ind(i,offset+X(i,s)) = 1;
    end
    wTrue(offset+1:offset+nStates(s),1) = (rand > .75)*randn(nStates(s),1);
    offset = offset+nStates(s);
end
y = X_ind*wTrue + randn(nInstances,1);

Xtrain = X_ind(1:floor(nInstances/2),:);
ytrain = y(1:floor(nInstances/2));
Xtest = X_ind(floor(nInstances/2)+1:end,:);
ytest = y(floor(nInstances/2)+1:end);

%% Set up variables and groups

% Initial guess of parameters
w = zeros(sum(nStates),1);
nVars = length(w);

% Set up groups
offset = 0;
groups = zeros(size(w));
for s = 1:length(nStates)
    groups(offset+1:offset+nStates(s),1) = s;
    offset = offset+nStates(s);
end
nGroups = max(groups);

%% Set up objective and compute maximum value of lambda

% Set up loss function
funObj1 = @(w)SquaredError(w,Xtrain,ytrain);

% Compute maximum value of lambda
[f,g] = funObj1(zeros(nVars,1));
grad_norms = sqrt(accumarray(groups(groups~=0),g(groups~=0).^2));
maxLambda = max(grad_norms);

% Set up validation function
testObj = @(w)SquaredError(w,Xtest,ytest);

%% Solve for regularization path

i = 1;
options.verbose = 0;
testNLL = testObj(w);
for lambda = maxLambda-maxLambda/100:-maxLambda/100:0
    
    % Make vector of regularization weights
    lambdaVect = lambda*ones(nGroups,1);
    
    % Set up regularization function
    funObj2 = @(w)groupL1regularizer(w,lambdaVect,groups);
    
    % Set up soft-threshold function
    funProj = @(w,stepSize)groupSoftThreshold(w,stepSize,lambdaVect,groups);
    
    i = i + 1;
    w(:,i) = minConf_BBST(funObj1,funObj2,w(:,i-1),funProj,options);
    testNLL(i) = testObj(w(:,i));
end
%% Show true parameters and parameters based on validation set
[minVal,minInd] = min(testNLL);
wMinTestNLL = w(:,minInd);
[wTrue wMinTestNLL]

%% Show full regularization path
imagesc(w);
colormap gray
title('Regularization Path');
w

