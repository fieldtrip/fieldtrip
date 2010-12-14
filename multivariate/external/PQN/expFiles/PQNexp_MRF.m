function [] = UGMlearn_Sachs(doBigData,nInstances,lambda,corrections,infer)

if doBigData
    load sachsDiscretizedData
else
    load 20news_w100
    [sorted sortInd]=sort(sum(documents'));
    data = 1+double(full(documents(sortInd(93:100),:)));
end

y = data';
y = y(1:nInstances,:);

[nInstances,nNodes] = size(y);
nFeatures = 0;
X = zeros(nInstances,nFeatures,nNodes);
adjInit = fixed_Full(nNodes);
nStates = max(y(:));

ising = 0;
tied = 0;
useMex = 1;

% Make edgeStruct
edgeStruct = UGM_makeEdgeStruct(adjInit,nStates,useMex);
nEdges = edgeStruct.nEdges;
Xedge = zeros(nInstances,nFeatures,nEdges);

% Add Biases
Xnode = [ones(nInstances,1,nNodes) X];
Xedge = [ones(nInstances,1,nEdges) Xedge];

% Make infoStruct and initialize Weights
infoStruct = UGM_makeInfoStruct(Xnode,Xedge,edgeStruct,ising,tied);
[w,v] = UGM_initWeights(infoStruct,@zeros);
nVars = numel(w)+numel(v);

% Make Groups
nodeGroups = zeros(size(w));
edgeGroups = zeros(size(v));
for e = 1:nEdges
    edgeGroups(:,:,e)=e;
end
groups = [nodeGroups(:);edgeGroups(:)];

% Make Objective Function
switch infer
    case 'pseudo'
        fprintf('Pseudolikelihood\n');
        funObj_sub = @(wv)UGM_PseudoLoss(wv,Xnode,Xedge,y,edgeStruct,infoStruct);
    case 'loopy'
        fprintf('Loopy\n');
        funObj_sub = @(wv)UGM_MRFLoss(wv,y,edgeStruct,infoStruct,@UGM_Infer_LBP);
    case 'exact'
        fprintf('Exact\n');
        funObj_sub = @(wv)UGM_MRFLoss(wv,y,edgeStruct,infoStruct,@UGM_Infer_Exact);
end

funObj = @(w)traceGroupL1(w,lambda,groups,funObj_sub);
global fValues;

%% Maximum Likelihood
if 0
    wv = minFunc(funObj,[w(:);v(:)]);
end

%% L1-Regularization
if 0
    lambdaNode = zeros(size(w));
    lambdaEdge = lambda*ones(size(v));
    options.order = -1;
    options.maxIter = 1000;
    wv = L1GeneralProjectedSubGradient(funObj,[w(:);v(:)],[lambdaNode(:);lambdaEdge(:)],options);
end

%% Group L1-Regularization

global fValues;
options.maxIter = 1000;

% Group-Grafting
fprintf('Running Grafting\n');
graftOptions.method = 'lbfgs';
graftOptions.maxIter = 1000;
fValues = [];
wvGraft = L1groupGraft(funObj,[w(:);v(:)],groups,lambda,graftOptions);
fGraft = fValues;

% SPG
fprintf('Running SPG\n');
options.normType = 2;
options.mode = 'spg';
fValues = [];
wvSPG = L1groupMinConF(funObj,[w(:);v(:)],groups,lambda,options);
fSPG = fValues;

% PQN
fprintf('Running PQN on Penalty Formulation\n');
options.mode = 'sop';
options.method = 'lbfgs';
options.corrections = corrections;
fValues = [];
wvPQN = L1groupMinConF(funObj,[w(:);v(:)],groups,lambda,options);
fPQN = fValues;

% Compute tau
if fPQN(end) < fSPG(end)
    [wtemp,vtemp] = UGM_splitWeights(wvPQN,infoStruct);
else
    [wtemp,vtemp] = UGM_splitWeights(wvSPG,infoStruct);
end
tau = sum(sqrt(accumarray(edgeGroups(:),vtemp(:).^2)));

% PQN (direct)
fprintf('Running PQN on Constrained Formulation\n');
funProj = @(w)groupL2Proj(w,tau,groups);
options.bbInit = 0;
options.SPGoptTol = 1e-6;
options.SPGiters = 500;
options.maxProject = inf;
options.SPGtestOpt = 1;
fValues = [];
wvPQN2 = minConF_PQN(funObj,[w(:);v(:)],funProj,options);
fPQN2 = fValues;

fileName = sprintf('PQNexp_MRF_%d_%d_%f_%d_%s.mat',doBigData,nInstances,lambda,corrections,infer);
save(fileName);
