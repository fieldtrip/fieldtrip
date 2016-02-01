function fit = glmnet(x, y, family, options)

%--------------------------------------------------------------------------
% glmnet.m: fit an elasticnet model path
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Fit a regularization path for the elasticnet at a grid of values for
%    the regularization parameter lambda. Can deal with all shapes of data.
%    Fits linear, logistic and multinomial regression models.
%
% USAGE: 
%    fit = glmnet(x, y)
%    fit = glmnet(x, y, family, options)
%
% EXTERNAL FUNCTIONS:
% options         = glmnetSet;                  provided with glmnet.m
%
% INPUT ARGUMENTS:
% x           Input matrix, of dimension nobs x nvars; each row is an
%             observation vector. Can be in sparse column format.
% y           Response variable. Quantitative for family =
%             'gaussian'. For family = 'binomial' should be either a vector
%             of two levels, or a two-column matrix of counts or
%             proportions. For family = 'multinomial', can be either a
%             vector of nc>=2 levels, or a matrix with nc columns of counts
%             or proportions.
% family      Reponse type. (See above). Default is 'gaussian'. 
% options     A structure that may be set and altered by glmnetSet (type 
%             help glmnetSet).
%
% OUTPUT ARGUMENTS:
% fit         A structure.
% fit.a0      Intercept sequence of length length(fit.lambda). 
% fit.beta    For "elnet" and "lognet" models, a nvars x length(lambda) 
%             matrix of coefficients. For "multnet", a list of nc such
%             matrices, one for each class.
% fit.lambda  The actual sequence of lambda values used.
% fit.dev     The fraction of (null) deviance explained (for "elnet", this 
%             is the R-square).
% fit.nulldev Null deviance (per observation).
% fit.df      The number of nonzero coefficients for each value of lambda.
%             For "multnet", this is the number of variables with a nonzero
%             coefficient for any class.
% fit.dfmat   For "multnet" only. A matrix consisting of the number of
%             nonzero coefficients per class.
% fit.dim     Dimension of coefficient matrix (ices).
% fit.npasses Total passes over the data summed over all lambda values.
% fit.jerr    Error flag, for warnings and errors (largely for internal 
%             debugging).
% fit.class   Type of regression - internal usage.
%
% DETAILS:
%    The sequence of models implied by lambda is fit by coordinate descent.
%    For family = 'gaussian' this is the lasso sequence if alpha = 1, else
%    it is the elasticnet sequence. For family = 'binomial' or family =
%    "multinomial", this is a lasso or elasticnet regularization path for
%    fitting the linear logistic or multinomial logistic regression paths.
%    Sometimes the sequence is truncated before options.nlambda values of
%    lambda have been used, because of instabilities in the logistic or
%    multinomial models near a saturated fit. glmnet(..., family =
%    'binomial') fits a traditional logistic regression model for the
%    log-odds. glmnet(..., family = 'multinomial') fits a symmetric
%    multinomial model, where each class is represented by a linear model
%    (on the log-scale). The penalties take care of redundancies. A
%    two-class "multinomial" model will produce the same fit as the
%    corresponding "binomial" model, except the pair of coefficient
%    matrices will be equal in magnitude and opposite in sign, and half the
%    "binomial" values. Note that the objective function for
%    "gaussian" is 
%                1 / (2 * nobs) RSS + lambda * penalty
%    , and for the logistic models it is 
%                1 / nobs - loglik + lambda * penalty
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
%    glmnetSet, glmnetPrint, glmnetPlot, glmnetPredict and glmnetCoef methods.
% 
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    g2=randsample(2,100,true);
%    g4=randsample(4,100,true);
%    fit1=glmnet(x,y);
%    glmnetPrint(fit1);
%    glmnetCoef(fit1,0.01) % extract coefficients at a single value of lambda
%    glmnetPredict(fit1,'response',x(1:10,:),[0.01,0.005]') % make predictions
%    fit2=glmnet(x,g2,'binomial');
%    fit3=glmnet(x,g4,'multinomial');
%
% DEVELOPMENT: 
%    14 Jul 2009: Original version of glmnet.m written.
%    26 Jan 2010: Fixed a bug in the description of y, pointed out by
%                 Peter Rijnbeek from Erasmus University.
%    09 Mar 2010: Fixed a bug of printing "ka = 2", pointed out by 
%                 Ramon Casanova from Wake Forest University.
%    25 Mar 2010: Fixed a bug when p > n in multinomial fitting, pointed
%                 out by Gerald Quon from University of Toronto

% Check input arguments
if nargin < 2
    error('more input arguments needed.');
end

if nargin < 3
    family = 'gaussian';
end

if nargin < 4
    options = glmnetSet;
end

% Prepare parameters
nlam = options.nlambda;

[nobs,nvars] = size(x);

weights = options.weights;
if isempty(weights)
    weights = ones(nobs,1);
end

maxit = options.maxit;

if strcmp(family, 'binomial') || strcmp(family, 'multinomial')
    [noo,nc] = size(y);
    kopt = double(options.HessianExact);
    if noo ~= nobs
        error('x and y have different number of rows');
    end
    if nc == 1
        classes = unique(y);
        nc = length(classes);
        indexes = eye(nc);
        y = indexes(y,:);
    end
    if strcmp(family, 'binomial')
        if nc > 2
            error ('More than two classes; use multinomial family instead');
        end
        nc = 1; % for calling multinet
    end
    if ~isempty(weights)
        % check if any are zero
        o = weights > 0;
        if ~all(o) %subset the data
            y = y(o,:);
            x = x(o,:);
            weights = weights(o);
            nobs = sum(o);
        end
        [my,ny] = size(y);
        y = y .* repmat(weights,1,ny);
    end
    % Compute the null deviance
    prior = sum(y,1);
    sumw = sum(sum(y));
    prior = prior / sumw;
    nulldev = -2 * sum(sum(y .* (ones(nobs, 1) * log(prior)))) / sumw;
elseif strcmp(family, 'gaussian')
    % Compute the null deviance
    ybar = y' * weights/ sum(weights);
    nulldev = (y' - ybar).^2 * weights / sum(weights);
    if strcmp(options.type, 'covariance')
        ka = 1;
    elseif strcmp(options.type, 'naive')
        ka = 2;
    else
        error('unrecognized type');
    end
else
    error('unrecognized family');
end

ne = options.dfmax;
if ne == 0
    ne = nvars + 1;
end
nx = options.pmax;
if nx == 0
    nx = min(ne * 1.2, nvars);
end
exclude = options.exclude;
if ~isempty(exclude)
    exclude = unique(exclude);
    if ~all(exclude > 0 & exclude <= nvars)
        error('Some excluded variables out of range');
    end
    jd = [length(exclude); exclude];
else
    jd = 0;
end
vp = options.penalty_factor;
if isempty(vp)
    vp = ones(nvars,1);
end
isd = double(options.standardize);
thresh = options.thresh;
lambda = options.lambda;
lambda_min = options.lambda_min;
if lambda_min == 0
    if nobs < nvars
        lambda_min = 5e-2;
    else
        lambda_min = 1e-4;
    end
end
if isempty(lambda)
    if (lambda_min >= 1)
        error('lambda_min should be less than 1');
    end
    flmin = lambda_min;
    ulam = 0;
else
    flmin = 1.0;
    if any(lambda < 0)
        error ('lambdas should be non-negative');
    end
    ulam = -sort(-lambda);
    nlam = length(lambda);
end

parm = options.alpha;

if strcmp(family, 'gaussian')
    [a0,ca,ia,nin,rsq,alm,nlp,jerr] = glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,weights,ka);
else
    [a0,ca,ia,nin,dev,alm,nlp,jerr] = glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,nc,maxit,kopt);
end

% Prepare output
lmu = length(alm);
ninmax = max(nin);
lam = alm;
if isempty(options.lambda)
    lam = fix_lam(lam); % first lambda is infinity; changed to entry point
end
errmsg = err(jerr, maxit, nx);
if errmsg.n == 1
    error(errmsg.msg);
elseif errmsg.n == -1
    warning(errmsg.msg);
end

if strcmp(family, 'multinomial')
    beta_list = {};
    a0 = a0 - repmat(mean(a0), nc, 1);
    dfmat=a0;
    dd=[nvars, lmu];
    if ninmax > 0
        ca = reshape(ca, nx, nc, lmu);
        ca = ca(1:ninmax,:,:);
        ja = ia(1:ninmax);
        [ja1,oja] = sort(ja);
        df = any(abs(ca) > 0, 2);
        df = sum(df, 1);
        df = df(:);
        for k=1:nc
            ca1 = reshape(ca(:,k,:), ninmax, lmu);
            cak = ca1(oja,:);
            dfmat(k,:) = sum(sum(abs(cak) > 0));
            beta = zeros(nvars, lmu);
            beta(ja1,:) = cak;
            beta_list{k} = beta;
        end
    else
        for k = 1:nc
            dfmat(k,:) = zeros(1,lmu);
            beta_list{k} = zeros(nvars, lmu);
        end
    end
    fit.a0 = a0;
    fit.beta = beta_list;
    fit.dev = dev;
    fit.nulldev = nulldev;
    fit.dfmat = dfmat;
    fit.df = df';
    fit.lambda = lam;
    fit.npasses = nlp;
    fit.jerr = jerr;
    fit.dim = dd;
    fit.class = 'multnet';
else
    dd=[nvars, lmu];
    if ninmax > 0
        ca = ca(1:ninmax,:);
        df = sum(abs(ca) > 0, 1);
        ja = ia(1:ninmax);
        [ja1,oja] = sort(ja);
        beta = zeros(nvars, lmu);
        beta (ja1, :) = ca(oja,:);
    else
        beta = zeros(nvars,lmu);
        df = zeros(1,lmu);
    end
    
    if strcmp(family, 'binomial')
        a0 = -a0;
        fit.a0 = a0;
        fit.beta = -beta; %sign flips make 2 arget class
        fit.dev = dev;
        fit.nulldev = nulldev;
        fit.df = df';
        fit.lambda = lam;
        fit.npasses = nlp;
        fit.jerr = jerr;
        fit.dim = dd;
        fit.class = 'lognet';
    else
        fit.a0 = a0;
        fit.beta = beta;
        fit.dev = rsq;
        fit.nulldev = nulldev;
        fit.df = df';
        fit.lambda = lam;
        fit.npasses = nlp;
        fit.jerr = jerr;
        fit.dim = dd;
        fit.class = 'elnet';
    end
end

%------------------------------------------------------------------
% End function glmnet
%------------------------------------------------------------------

function new_lam = fix_lam(lam)

new_lam = lam;
llam=log(lam);
new_lam(1)=exp(2*llam(2)-llam(3));

%------------------------------------------------------------------
% End private function fix_lam
%------------------------------------------------------------------

function output = err(n,maxit,pmax)

if n==0
    output.n=0;
    output.msg='';
elseif n>0 %fatal error
    if n<7777
        msg='Memory allocation error; contact package maintainer';
    elseif n==7777 
        msg='All used predictors have zero variance';
    elseif (8000<n) && (n<9000)
        msg=sprintf('Null probability for class %d < 1.0e-5', n-8000);
    elseif (9000<n) && (n<10000)
        msg=sprintf('Null probability for class %d > 1.0 - 1.0e-5', n-9000);
    elseif n==10000
        msg='All penalty factors are <= 0';
    end
    output.n=1
    output.msg=['in glmnet fortran code - %s',msg];
elseif n<0 %non fatal error
    if n > -10000 
        msg=sprintf('Convergence for %dth lambda value not reached after maxit=%d iterations; solutions for larger lambdas returned', -n, maxit);
    elseif n < -10000 
        msg=sprintf('Number of nonzero coefficients along the path exceeds pmax=%d at %dth lambda value; solutions for larger lambdas returned', pmax, -n-10000);
    end
    output.n=-1;
    output.msg=['from glmnet fortran code - ',msg];
end

%------------------------------------------------------------------
% End private function err
%------------------------------------------------------------------
