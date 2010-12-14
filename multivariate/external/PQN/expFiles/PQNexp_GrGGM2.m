function [] = PQNexp_GrGGM2(doBigData,nfcv)

if doBigData
    load genes.mat
    load grouping.mat
    tikhonovStrength = .2;
else
    load 20news_w100
    X = double(documents);
    assignWords
    tikhonovStrength = 1e-4;
end

% Get groups and compute number of elements per group
groups = z;
[z,idx] = sort(groups);
z(2:end) = z(2:end) - z(1:end-1);
z = find(z);
z = z(2:end) - z(1:end-1);
z(end+1) = length(groups) - sum(z);
nIndices = makeIndices(groups); % used by projection

% Reorder the data to have consecutive groups
groups = groups(idx);
X      = X(idx,:);
X      = X - repmat(mean(X,2),1,size(X,2));
X      = spdiags(1./std(X,0,2),0,size(X,1),size(X,1)) * X;

% Parameters of optimizers
options.method = 'lbfgs';
options.corrections = 5;

[p,n] = size(X);

% =====================================================================
% L1 and Block-L1 code
% =====================================================================

lambdaVec = logspace(-4, -1.5, 50);

d = ceil(size(X,2) / nfcv);
resultL1   = zeros(nfcv,length(lambdaVec));
resultLInf = zeros(nfcv,length(lambdaVec));
resultL2   = zeros(nfcv,length(lambdaVec));

options.verbose = 0;
for i=1:nfcv;
    fprintf('Cross-validation fold %d of %d\n',i,nfcv)
    
    % Get indices
    perm = randperm(n);
    idx1 = perm(1:floor(4*n/5));
    idx2 = perm(floor(4*n/5)+1:end);

    Xbar   = X(:,idx1);
    Sigma1 = Xbar*Xbar' / size(Xbar,2);
    Sigma1 = Sigma1 + tikhonovStrength * eye(size(Sigma1,1));
    Xbar   = X(:,idx2);
    Sigma2 = Xbar*Xbar' / size(Xbar,2);
    
    funObj = @(X)logdetFunction(X,Sigma1);

    K = inv(Sigma1);
    resultTik(i) = -667*log(2*pi)/2 + logdet(K,-Inf) / 2 - traceMatProd(Sigma2,K) / 2;
    
    %% L1
    W = [];
    for j=1:length(lambdaVec)
        lambdaBar = lambdaVec(j);
        lambdaFull = ones(p)*lambdaBar;

        fprintf('L1 with lambda = %f\n',lambdaBar);
        if 1
            X_init = diag(diag(lambdaFull));
            funProj = @(X)boundProject(X,-lambdaFull(:),lambdaFull(:));
            W = minConF_SPG(funObj,X_init(:),funProj,options);
            K = inv(Sigma1+reshape(W,[p p]));
        else
            [K,W] = Algorithm1(Sigma1,lambdaFull,W,Sigma2);
        end
        resultL1(i,j) = -667*log(2*pi)/2 + logdet(K,-Inf) / 2 - traceMatProd(Sigma2,K) / 2;

        save(sprintf('PQNexp_GrGGM2_%d_%d.mat',doBigData,nfcv),...
            'resultTik','resultL1','resultL2','resultLInf','lambdaVec');
    end

    %% L1-inf
    W = [];
    for j=1:length(lambdaVec)

        lambdaBar = lambdaVec(j);

        % Lambda is scalar times number of elements per group
        lambda = lambdaBar * z' * z;
        lambda = lambda - diag(diag(lambda));
        lambda = lambda + diag(ones(length(lambda),1)*lambdaBar);
        lambda = lambda(:);

            fprintf('Group L{1,inf} with lambda = %f\n',lambdaBar);
        if 1
            X_init = lambdaBar*eye(p);
            funProj = @(W) projectLinf1(W,p,nIndices,lambda);
            W = minConF_SPG(funObj,X_init(:),funProj,options);
            K = inv(Sigma1+reshape(W,[p p]));
        else
            [K,W,f] = Algorithm3BlockMatrix(Sigma1, groups, lambda, Inf, W, Sigma2);
        end

        resultLInf(i,j) = -667*log(2*pi)/2 + logdet(K,-Inf) / 2 - traceMatProd(Sigma2,K) / 2;

                save(sprintf('PQNexp_GrGGM2_%d_%d.mat',doBigData,nfcv),...
            'resultTik','resultL1','resultL2','resultLInf','lambdaVec');
    end

    %% L1-2
    W = [];
    for j=1:length(lambdaVec)
        lambdaBar = lambdaVec(j);

        % Lambda is scalar times number of elements per group
        lambda = lambdaBar * z' * z;
        lambda = lambda - diag(diag(lambda));
        lambda = lambda + diag(ones(length(lambda),1)*lambdaBar);
        lambda = lambda(:);

        fprintf('Group L{1,2} with lambda = %f\n',lambdaBar);
        if 1
            X_init = lambdaBar*eye(p);
        funProj = @(W) projectLinf2(W,p,nIndices,lambda);
        W = minConF_SPG(funObj,X_init(:),funProj,options);
        K = inv(Sigma1+reshape(W,[p p]));
        else
        [K,W,f] = Algorithm3BlockMatrix(Sigma1, groups, lambda, 2, W, Sigma2);
        end
        
        
        resultL2(i,j) = -667*log(2*pi)/2 + logdet(K,-Inf) / 2 - traceMatProd(Sigma2,K) / 2;

                save(sprintf('PQNexp_GrGGM2_%d_%d.mat',doBigData,nfcv),...
            'resultTik','resultL1','resultL2','resultLInf','lambdaVec');
    end

end % For i


