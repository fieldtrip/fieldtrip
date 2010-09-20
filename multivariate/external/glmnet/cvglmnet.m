function CVerr = cvglmnet(x,y,nfolds,foldid,type,family,options,verbous)
% Do crossvalidation of glmnet model. The coordinate descent algorithm 
% chooses a set of lambda to optimize over based on the input given in
% the options struct.Parameter tuning should be done in such a way that for
% a fixed alpha, a set of lambda values are evaluated. Basically it does
% not matter if lambda corresponds across alpha values, as each cv-plot
% should inspected seperatly.
% So e.g. to find optimal tuning parameters, fit a total of 10 alpha
% values, beeing alpha = 0:0.1:1.lambdas will then be chosen according to 
% the specific alpha. 
% Call: CVerr = cvglmnet(x,y,nfolds,foldid,type,family,options,verbous)
% Example:
% x=randn(100,2000);
% y=randn(100,1);
% g2=randsample(2,100,true);
% CVerr=cvglmnet(x,y,100,[],'response','gaussian',glmnetSet,1);
% CVerr=cvglmnet(x,g2,100,[],'response','binomial',glmnetSet,1);
% x         : Covariates
% y         : Response (For now only elastic net = continous data supported
% nfolds    : How many folds to evaluate. nfolds = size(x,1) for LOOCV
% foldid    : Possibility for supplying own folding series. [] for nothing
% type      : Used in the glmnetPredict routine. (Now only "response" works)
% family    : Used in glmnet routine. (Now only "gaussian" works)
% options   : See function glmnetSet()
% verbous   : Print model plot
% 
% Written by Bjørn Skovlund Dissing (27-02-2010)
glmnet_object = glmnet(x, y, family,options);
options.lambda = glmnet_object.lambda;
options.nlambda = length(options.lambda);
N = size(x,1);
if (isempty(foldid))
    foldid = randsample([repmat(1:nfolds,1,floor(N/nfolds)) 1:mod(N,nfolds)],N);
else
    nfolds = max(foldid)
end
predmat = glmnetPredict(glmnet_object, type,x, options.lambda);

for i=1:nfolds
    which=foldid==i;
    if verbous, disp(['Fitting fold # ' num2str(i) ' of ' num2str(nfolds)]);end
    cvfit = glmnet(x(~which,:), y(~which),family, options);
    predmat(which,:) = glmnetPredict(cvfit, type,x(which,:),options.lambda);
end
yy=repmat(y,1,length(options.lambda));
if strcmp(family,'gaussian')
    cvraw=(yy-predmat).^2;
elseif strcmp(family,'binomial')
    if     strcmp(type,'response')
        cvraw=-2*((yy==2).*log(predmat)+(yy==1).*log(1-predmat));
    elseif strcmp(type,'class')
        cvraw=double(yy~=predmat);
    end
elseif strcmp(family,'multinomial')
    error('Not implemented yet')
end
CVerr.cvm=mean(cvraw);
CVerr.stderr=sqrt(var(cvraw)/N);
CVerr.cvlo=CVerr.cvm-CVerr.stderr;
CVerr.cvup=CVerr.cvm+CVerr.stderr;
% if there are several minima, choose largest lambda of the smallest cvm
CVerr.lambda_min=max(options.lambda(CVerr.cvm<=min(CVerr.cvm)));
%Find stderr for lambda(min(sterr))
semin=CVerr.cvup(options.lambda==CVerr.lambda_min);
% find largest lambda which has a smaller mse than the stderr belonging to
% the largest of the lambda belonging to the smallest mse
% In other words, this defines the uncertainty of the min-cv, and the min
% cv-err could in essence be any model in this interval.
CVerr.lambda_1se=max(options.lambda(CVerr.cvm<semin));
CVerr.glmnetOptions=options;
CVerr.glmnet_object = glmnet_object;
if verbous, cvglmnetPlot(CVerr);end
end