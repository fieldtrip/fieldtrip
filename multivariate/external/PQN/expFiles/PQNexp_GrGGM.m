function [] = L1precisionGroup_PQN(doBigData,lambdaBar,corrections)

fileName = sprintf('PQNexp_GrGGM_%d_%f_%d.mat',doBigData,lambdaBar,corrections);

%% Data Set

if doBigData
    load genes.mat
    load grouping.mat
else
    load 20news_w100
    X = double(documents);
    assignWords
end

global trace
trace = 1;
global fValues

%% Preprocessing for groups

[p,n] = size(X);

% Get groups and compute number of elements per group
groups = z;
[z,idx] = sort(groups);
z(2:end) = z(2:end) - z(1:end-1);
z = find(z);
z = z(2:end) - z(1:end-1);
z(end+1) = length(groups) - sum(z);

% Reorder the data to have consecutive groups
groups = groups(idx);
X      = X(idx,:);
XstdCol = standardizeCols(X')';
X      = X - repmat(mean(X,2),1,size(X,2));
X      = spdiags(1./std(X,0,2),0,size(X,1),size(X,1)) * X;

Xbar = X;
Sigma1 = Xbar*Xbar' / size(Xbar,2);

% Lambda is scalar times number of elements per group
lambda = lambdaBar * z' * z;
lambda = lambda - diag(diag(lambda));
lambda = lambda + diag(ones(length(lambda),1)*lambdaBar);
lambda = lambda(:);

fprintf('Running PG\n');
fValues = [];
wPG = Algorithm3BlockMatrix(Sigma1, groups, lambda, Inf);
fPG = fValues;
save(fileName);

%% Setup L{inf,1} Projection
groups  = groups(:);
nGroups = max(groups);
indices = cell(nGroups,1);
for i=1:nGroups
    indices{i} = find(groups == i);
end
nIndices = zeros(size(indices));
for i=1:nGroups
    nIndices(i) = length(indices{i});
end
funProj = @(W) projectLinf1(W,p,nIndices,lambda); % (W, nIndices, lambda)
%%

X_init = lambdaBar*eye(p);
funObj = @(X)logdetFunction(X,Sigma1);
options = [];
fValues = [];
fprintf('Running SPG\n');
W = minConF_SPG(funObj,X_init(:),funProj,options);
fSPG = fValues;
wSPG = inv(Sigma1+reshape(W,[p p]));
save(fileName);

options.bbInit = 0;
options.SPGoptTol = 1e-16;
options.SPGiters = 500;
options.maxProject = inf;
options.SPGtestOpt = 1;
options.corrections = corrections;
fValues = [];
W = minConF_PQN(funObj,X_init(:),funProj,options);
fPQN = fValues;
wPQN = inv(Sigma1+reshape(W,[p p]));
save(fileName);