function [] = L1precision_PQN(doBigData,lambda,corrections)

fileName = sprintf('PQNexp_GGM_%d_%f_%d.mat',doBigData,lambda,corrections);

%% Data Set

if doBigData
    load genes.mat
    X = standardizeCols(X');
    [n,p] = size(X);
    sigma = (1/n)*(X'*X);
else
    load 20news_w100
    X = double(documents)';
    [n,p] = size(X);
    X = standardizeCols(X);
    sigma = (1/n)*(X'*X);
end

global trace
trace = 1;
global fValues

%%

fprintf('Running BCD\n');
fValues = [];
wBCD = L1precisionBCD_traced(sigma,lambda);
fBCD = fValues;
save(fileName);

lambdaFull = lambda*ones(p);

fprintf('Running PG\n');
fValues = [];
wPG = Algorithm1(sigma,lambdaFull);
fPG = fValues;
save(fileName);

fprintf('Running SPG\n');
X_init = diag(diag(lambdaFull));
funObj = @(X)logdetFunction(X,sigma);
funProj = @(X)boundProject(X,-lambdaFull(:),lambdaFull(:));
options = [];
fValues = [];
W = minConF_SPG(funObj,X_init(:),funProj,options);
fSPG = fValues;
wSPG = inv(sigma+reshape(W,[p p]));
save(fileName);

fprintf('Running PQN\n');
options.bbInit = 0;
options.SPGoptTol = 1e-16;
options.SPGiters = 500;
options.maxProject = inf;
options.SPGtestOpt = 1;
options.corrections = corrections;
fValues = [];
W = minConF_PQN(funObj,X_init(:),funProj,options);
fPQN = fValues;
wPQN = inv(sigma+reshape(W,[p p]));
save(fileName);
