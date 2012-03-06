function result = glmnetPredict(object, type, newx, s)

%--------------------------------------------------------------------------
% glmnetPredict.m: make predictions from a "glmnet" object.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Similar to other predict methods, this functions predicts fitted
%    values, logits, coefficients and more from a fitted "glmnet" object.
%
% USAGE: 
%    glmnetPredict(object)
%    glmnetPredict(object, type)
%    glmnetPredict(object, type, newx)
%    glmnetPredict(object, type, newx, s)
%
% INPUT ARGUMENTS:
% fit         Fitted "glmnet" model object.
% type        Type of prediction required. Type "link" gives the linear
%             predictors for "binomial" or "multinomial" models; for
%             "gaussian" models it gives the fitted values. Type "response"
%             gives the fitted probabilities for "binomial" or
%             "multinomial"; for "gaussian" type "response" is equivalent
%             to type "link". Type "coefficients" computes the coefficients
%             at the requested values for s. Note that for "binomial"
%             models, results are returned only for the class corresponding
%             to the second level of the factor response. Type "class"
%             applies only to "binomial" or "multinomial" models, and
%             produces the class label corresponding to the maximum
%             probability. Type "nonzero" returns a list of the indices of
%             the nonzero coefficients for each value of s.
% newx        Matrix of new values for x at which predictions are to be 
%             made. Must be a matrix; This argument is not used for 
%             type=c("coefficients","nonzero")
% s           Value(s) of the penalty parameter lambda at which predictions
%             are required. Default is the entire sequence used to create
%             the model.
%
% DETAILS:
%    The shape of the objects returned are different for "multinomial"
%    objects. glmnetCoef(fit, ...) is equivalent to glmnetPredict(fit, "coefficients", ...)
%
% LICENSE: GPL-2
%
% DATE: 14 Jul 2009
%
% AUTHORS:
%    Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
%    Fortran code was written by Jerome Friedman 
%    R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
%    MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
%    Department of Statistics, Stanford University, Stanford, California, USA.
%
% REFERENCES:
%    Friedman, J., Hastie, T. and Tibshirani, R. (2009)
%    Regularization Paths for Generalized Linear Models via Coordinate Descent.
%    Journal of Statistical Software, 33(1), 2010
%
% SEE ALSO:
%    glmnet, glmnetSet, glmnetPrint, glmnetPlot and glmnetCoef methods.
% 
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    g2=randsample(2,100,true);
%    g4=randsample(4,100,true);
%    fit1=glmnet(x,y);
%    glmnetPredict(fit1,'link',x(1:5,:),[0.01,0.005]') % make predictions
%    glmnetPredict(fit1,'coefficients')
%    fit2=glmnet(x,g2,'binomial');
%    glmnetPredict(fit2, 'response', x(2:5,:))
%    glmnetPredict(fit2, 'nonzero')
%    fit3=glmnet(x,g4,'multinomial');
%    glmnetPredict(fit3, 'response', x(1:3,:), 0.01)
%
% DEVELOPMENT: 
%    14 Jul 2009: Original version of glmnet.m written.
%    20 Oct 2009: Fixed a bug in bionomial response, pointed out by Ramon
%                 Casanov from Wake Forest University.
%    26 Jan 2010: Fixed a bug in multinomial link and class, pointed out by
%                 Peter Rijnbeek from Erasmus University.

if nargin < 2
    type = 'link';
end

if nargin < 3
    newx = [];
end

if nargin < 4
    s = object.lambda;
end

if strcmp(object.class, 'elnet')
    a0=transpose(object.a0);
    nbeta=[a0; object.beta];
    if nargin == 4
        lambda=object.lambda;
        lamlist=lambda_interp(lambda,s);
        nbeta=nbeta(:,lamlist.left).*repmat(lamlist.frac',size(nbeta,1),1) +nbeta(:,lamlist.right).*(1-repmat(lamlist.frac',size(nbeta,1),1));
    end
    if strcmp(type, 'coefficients')
        result = nbeta;
    elseif strcmp(type, 'link')
        result = [ones(size(newx,1),1), newx] * nbeta;
    elseif strcmp(type, 'response')
        result = [ones(size(newx,1),1), newx] * nbeta;
    elseif strcmp(type, 'nonzero')
        result = nonzeroCoef(nbeta(2:size(nbeta,1),:), true);
    else
        error('Unrecognized type');
    end
elseif strcmp(object.class, 'lognet')
    
    a0=transpose(object.a0);
    nbeta=[object.a0; object.beta];
    if nargin == 4    
        lambda=object.lambda;
        lamlist=lambda_interp(lambda,s);
        nbeta=nbeta(:,lamlist.left).*repmat(lamlist.frac',size(nbeta,1),1) +nbeta(:,lamlist.right).*(1-repmat(lamlist.frac',size(nbeta,1),1));
    end
    %%%   remember that although the fortran lognet makes predictions
    %%%   for the first class, we make predictions for the second class
    %%%   to avoid confusion with 0/1 responses.
    %%%   glmnet flipped the signs of the coefficients 
    if strcmp(type,'coefficients')
        result = nbeta;
    elseif strcmp(type,'nonzero')
        result = nonzeroCoef(nbeta(2:size(nbeta,1),:), true);
    else
        nfit = [ones(size(newx,1),1), newx] * nbeta;
        
        if strcmp(type,'response')
            pp=exp(-nfit);
            result = 1./(1+pp);
        elseif strcmp(type,'link')
            result = nfit;
        elseif strcmp(type,'class')
            result = (nfit > 0) * 2 + (nfit <= 0) * 1;
        else
            error('Unrecognized type');
        end
    end
elseif strcmp(object.class, 'multnet')
    a0=object.a0;
    nbeta=object.beta;
    nclass=size(a0,1);
    nlambda=length(s);
    if nargin == 4
        lambda=object.lambda;
        lamlist=lambda_interp(lambda,s);
        for i=1:nclass
            kbeta=[a0(i,:); nbeta{i}];
            kbeta=kbeta(:,lamlist.left)*lamlist.frac +kbeta(:,lamlist.right)*(1-lamlist.frac);
            nbeta{i}=kbeta;
        end
    else
        for i=1:nclass
            nbeta{i} = [a0(i,:);nbeta{i}];
        end
    end
    if strcmp(type, 'coefficients')
        result = nbeta;
    elseif strcmp(type, 'nonzero')
        for i=1:nclass
            result{i}=nonzeroCoef(nbeta{i}(2:size(nbeta{i},1),:),true);
        end
    else    
        npred=size(newx,1);
        dp = zeros(nclass,nlambda,npred);
        for i=1:nclass
            fitk = [ones(size(newx,1),1), newx] * nbeta{i};
            dp(i,:,:)=dp(i,:,:)+reshape(transpose(fitk),1,nlambda,npred);
        end
        if strcmp(type, 'response')
            pp=exp(dp);
            psum=sum(pp,1);
            result = permute(pp./repmat(psum,nclass,1),[3,1,2]);
        elseif strcmp(type, 'link')
            result=permute(dp,[3,1,2]);
        elseif strcmp(type, 'class')
            dp=permute(dp,[3,1,2]);
            result = [];
            for i=1:size(dp,3)
                result = [result, softmax(dp(:,:,i))];
            end
        else
            error('Unrecognized type');
        end
    end
else
    error('Unrecognized class');
end

%-------------------------------------------------------------
% End private function glmnetPredict
%-------------------------------------------------------------

function result = lambda_interp(lambda,s)
% lambda is the index sequence that is produced by the model
% s is the new vector at which evaluations are required.
% the value is a vector of left and right indices, and a vector of fractions.
% the new values are interpolated bewteen the two using the fraction
% Note: lambda decreases. you take:
% sfrac*left+(1-sfrac*right)
  
  if length(lambda)==1 % degenerate case of only one lambda
    nums=length(s);
    left=ones(nums,1);
    right=left;
    sfrac=ones(nums,1);
else
      s(s > max(lambda)) = max(lambda);
      s(s < min(lambda)) = min(lambda);
      k=length(lambda);
      sfrac =(lambda(1)-s)/(lambda(1) - lambda(k));
      lambda = (lambda(1) - lambda)/(lambda(1) - lambda(k));
      coord = interp1(lambda, 1:length(lambda), sfrac);
      left = floor(coord);
      right = ceil(coord);
      sfrac=(sfrac-lambda(right))./(lambda(left) - lambda(right));
      sfrac(left==right)=1;
  end
  result.left = left;
  result.right = right;
  result.frac = sfrac;

%-------------------------------------------------------------
% End private function lambda_interp
%-------------------------------------------------------------
  
function result = softmax(x, gap) 
if nargin < 2
    gap = false;
end
d = size(x);
maxdist = x(:, 1);
pclass = repmat(1, d(1), 1);
for i =2:d(2)
    l = x(:, i) > maxdist;
    pclass(l) = i;
    maxdist(l) = x(l, i);
end
if gap
    x = abs(maxdist - x);
    x(1:d(1), pclass) = x * repmat(1, d(2));
    gaps = pmin(x);
end
if gap
    result = {pclass, gaps};
else
    result = pclass;
end

%-------------------------------------------------------------
% End private function softmax
%-------------------------------------------------------------
  