function [] = PQNexp_CRF(doBigData,lambda,nInstances,nVars)

fileName = sprintf('PQNexp_CRF_%d_%d_%d_%d.mat',doBigData,lambda,nInstances,nVars);

if doBigData
    load coNLL_train.mat
else
    load wordData.mat
end

% Make Problem smaller
X = X(:,1:nVars);
trainNdx = 1:nInstances;

% Set-up crfChain data structures
nFeatures = max(X);
nStates = max(y);
[w,v_start,v_end,v] = crfChain_initWeights(nFeatures,nStates,'zeros');
featureStart = cumsum([1 nFeatures(1:end)]);
sentences = crfChain_initSentences(y);
nSentences = size(sentences,1);
maxSentenceLength = 1+max(sentences(:,2)-sentences(:,1));

% Initialize parameter vector, set-up small L2-regularizer
w_init = [w(:);v_start(:);v_end(:);v(:)];
nVars = length(w_init);
lambda = lambda*ones(size(w_init));
lambdaL2 = 1e-4;

% Set up objective function
funObj_sub = @(wv)crfChain_lossL2C(wv,X,y,nStates,nFeatures,featureStart,sentences(trainNdx,:),maxSentenceLength,lambdaL2);

% Set up objective function wrapper that will keep track of objective
% values
funObj = @(w)traceL1(w,lambda,funObj_sub);
global fValues;

% Parameters of methods
options.order = -1;
options.optTol = 1e-10;
options.maxIter = 500;
options.corrections = 5;

fprintf('Bound-Constrained L-BFGS\n');
fValues = [];
wPr = L1GeneralProjection(funObj,w_init,lambda,options);
fPr = fValues;
save(fileName,'fPr');

fprintf('Orthant-Wise Descent\n');
fValues = [];
wOW = L1GeneralOrthantWise(funObj,w_init,lambda,options);
fOW = fValues;
save(fileName,'fPr','fOW');

% Compute value of tau that is equivalent to this value of lambda
if fPr(end) <= fOW(end)
    fprintf('Projection Found a Better Solution\n');
    tau = sum(abs(wPr));
else
    fprintf('Orthant-Wise Found a Better Solution\n');
    tau = sum(abs(wOW));
end

fprintf('SPG (constrained formulation, direct projeciton)\n');
fValues = [];
funProj = @(w)sign(w).*projectRandom2C(abs(w),tau);
wSPG = minConF_SPG(funObj,w_init,funProj,options);
fSPG = fValues;
save(fileName,'fPr','fOW','fSPG');

fprintf('Projected Quasi-Newton (constrained formulation, direct projeciton)\n');
options.bbInit = 0;
options.SPGoptTol = options.optTol;
options.SPGiters = options.maxIter;
options.maxProject = inf;
options.SPGtestOpt = 1;
fValues = [];
wPQN = minConF_PQN(funObj,w_init,funProj,options);
fPQN = fValues;
save(fileName,'fPr','fOW','fSPG','fPQN');

fprintf('Projected Quasi-Newton (penalty formulation, non-negative variables)\n');
fValues = [];
w_init = [w_init.*(w_init >= 0);-w_init.*(w_init <= 0)];
funObjWrap = @(w)nonNegGrad(w,lambda,funObj);
funProj = @(w)boundProject(w,zeros(2*nVars,1));
wPQNB = minConF_PQN(funObjWrap,w_init,funProj,options);
wPQNB = wPQNB(1:nVars)-wPQNB(nVars+1:end);
fPQNB = fValues;
save(fileName,'fPr','fOW','fSPG','fPQN','fPQNB');
