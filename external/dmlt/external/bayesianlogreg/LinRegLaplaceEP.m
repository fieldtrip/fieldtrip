function [Gauss,terms,logp,EP_ERR,logphist] = LinRegLaplaceEP(y,X,K,sig2,varargin)

% specialized version of our fast variant of
% Expectation Propagation for Logistic Regression (2 classes) with a (correlated) Laplace prior
% in the degenerate case, i.e., when the number of features >> number of samples
% (also works when number of samples > number of features, but then less stable!!!)
%
% input:
% labels    = an N x 1 vector of class labels [1,2]
% examples  = an N x M matrix of input data
% K         = the prior precision matrix of size M x M
%
% Note: the bias term should be explicitly added to K and the examples!
%
% regression parameters refers to betas
% auxiliary variables refers to u and v whose precision matrix auxK = inv(Lambda) couples the features.
%
% we have a precision matrix of the form
%
% | K_beta         |
% |        K_u     |
% |            K_v |
%
% priorGauss: struct with fields
%
%     hatK       diagonal of precision matrix (number of samples x 1);
%                   (initially zero)
%     diagK      diagonal of precision matrix (number of features x 1)
%                   (initially zero)
%     
%     precision matrix of regression parameters K_beta = A' hatK A + diagK
%
%     h          canonical mean (number of features x 1)
%                   (initially zero)
%
%     auxK       precision matrix of auxiliary variables (number of features x number of features; sparse)
%                   (contains the covariance structure of interest)
%
%     A          feature matrix (number of samples x number of features)
%
% terms: struct with fields
%
%     hatK       as priorGauss
%     diagK      as priorGauss
%     hath       canonical mean (number of samples x 1)
%     h          canonical mean (number of features x 1)
%                                canonical mean of regression parameters = h + A' hath
%     auxK       as priorGauss, but then only diagonal elements (number of features x 1)
%
% opt: struct with fields (all optional, defaults in brackets)
%
%     maxstepsize   maximum step size [1]
%     fraction      fraction or power for fractional/power EP [1]
%     niter         maximum number of iterations [100]
%     tol           convergence criterion [1e-5]
%     nweights      number of points for numerical integration [20]
%     temperature   temperature for simulated annealing [1]
%
% Gauss: struct with fields as in priorGauss plus
%
%     hatB       diagonal of projected covariance (number of samples x 1)
%                                                               hatB = A Covariance A'
%     hatn       projected mean (number of samples x 1)
%                                                               hatn = A m
%     diagC      diagonal of covariance matrix of regression parameters (number of features x 1)
%     auxC       diagonal of covariance matrix of auxiliary variables (number of features x 1)
%
%
% logp: estimated marginal loglikelihood (implementation doubtful...)
%
% comptime: computation time

fprintf('starting EP\n');
EP_ERR = 0;

tic

%% initialization

% parse opt  
%disp(varargin{:})
%length(varargin{:})

opt = [];
for i=1:2:(length(varargin)-1)
  opt.(varargin{i}) = varargin{i+1};
end

if ~isfield(opt,'maxstepsize'), opt.maxstepsize = 1; end
if ~isfield(opt,'fraction'),    opt.fraction = 0.95; end
if ~isfield(opt,'niter'),       opt.niter = 1000; end
if ~isfield(opt,'tol'),         opt.tol = 1e-5; end
if ~isfield(opt,'nweights'),    opt.nweights = 50; end
if ~isfield(opt,'temperature'), opt.temperature = 1; end
if ~isfield(opt,'lambda'),      opt.lambda = 1; end

stepsize = opt.maxstepsize;

[nobs,nfeatures] = size(X);

%% create priorGauss and terms


% construct Gaussian representation
v                = 1./sig2; 
priorGauss.A     = X;
priorGauss.hatK  = v*ones(nobs,1);
priorGauss.h     = v*X'*y;
priorGauss.diagK = zeros(nfeatures,1);
priorGauss.auxK  = K;

% compute the additional terms coming form the prior for the evidence
cholK  = chol(priorGauss.auxK,'lower');
LogPriorRestAux2 = (2*sum(log(full(diag(cholK))))-nfeatures*log(2*pi)); % no half because we need it two times
LogPriorRestLL   = 0.5*(-v*y'*y - nobs*log(2*pi) + nobs*log(v) ); 

% construct term representation

%terms.hatK  = zeros(nobs,1);
%terms.hath  = zeros(nobs,1);
terms.diagK = ones(nfeatures,1)*opt.lambda/10;
terms.auxK  = zeros(nfeatures,1);
terms.h     = zeros(nfeatures,1);

%% precompute points and weights for numerical integration

[xhermite,whermite] = gausshermite(opt.nweights);
xhermite = xhermite(:);       % nhermite x 1
whermite = whermite(:);       % nhermite x 1

[xlaguerre,wlaguerre] = gausslaguerre(opt.nweights);
xlaguerre = xlaguerre(:);     % nlaguerre x 1
wlaguerre = wlaguerre(:);     % nlaguerre x 1

% divide all canonical parameters by the temperature
if opt.temperature ~= 1,
   [priorGauss,terms] = correct_for_temperature(priorGauss,terms,opt.temperature);
end

%% build initial Gauss
Gauss        = update_Gauss(priorGauss,terms);               
Gauss        = canonical_to_moments(Gauss);
[myGauss,ok] = project_all(Gauss,terms,opt.fraction);
if ~ok,
   error('improper cavity distributions\n');
end

%% enter the iterations 

logp     = 0;
logpold  = 2*opt.tol;
change   = 0;
teller   = 0;
logphist = zeros(1,opt.niter);

while abs(logp-logpold) > stepsize*opt.tol && teller < opt.niter,
  
   teller    = teller+1;
   logpold   = logp;
   oldchange = change;
   
   % compute the new term approximation by applying an adf update on all the cavity approximations
   [fullterms,crosslogz] = ...
      adfupdate_all(myGauss,terms,xhermite,whermite,xlaguerre,wlaguerre,opt.fraction,opt.temperature);

   ok  = 0;
   ok1 = 1;
   ok2 = 1;
   
   while ~ok,
      % try to replace the old term approximations by the new term approximations and check whether the new Gauss is still fine
      
      % !!! THIS IS THE ORIGINAL LINE !!!
      
      % [newGauss,newterms,ok1] = try_update(Gauss,fullterms,terms,stepsize);
      
      % you might use this instead - sometimes it works better
      
      [newGauss,newterms,ok1] = try_update(Gauss,fullterms,terms,0.1);
     
      if ok1,
         % compute all cavity approximations needed for the next EP updates and check whether they are all fine         
         [newGauss,myGauss,logZappx,ok2] = try_project(newGauss,newterms,opt.fraction);
      end
      
      ok = (ok1 & ok2);
      
      if ok,    % accept
      
        terms    = newterms;
        Gauss    = newGauss;
        stepsize = min(opt.maxstepsize,stepsize*1.9);
         
      else  % try with smaller stepsize
      
        stepsize = stepsize/2;
         if ~ok1,
            fprintf('improper full covariance: lowering stepsize to %g\n',stepsize');
         elseif ~ok2,
            fprintf('improper cavity covariance: lowering stepsize to %g\n',stepsize');
         end
         if stepsize < 1e-10,
            warning('Cannot find an update that leads to proper cavity approximations.');
            EP_ERR = 1;
            teller = opt.niter;
            break;
         end

      end

   end
   
   % compute marginal moments
   
   if ok,
   
      % compute marginal loglikelihood
   
      %logp = sum(crosslogz) + logdet/2 + rest
      
      CorrTerm1 = (myGauss.m.^2)./myGauss.diagC + log(full(myGauss.diagC)) - (newGauss.m.^2)./newGauss.diagC - log(full(newGauss.diagC));
      CorrTerm2 = log(myGauss.auxC) - log(newGauss.auxC);
      CorrTerm  = 0.5*CorrTerm1 + 0.5*(CorrTerm2 + CorrTerm2); 
      logp      = sum(crosslogz + CorrTerm)./opt.fraction + logZappx + (LogPriorRestLL + LogPriorRestAux2);   

      logphist(teller) = logp;	
       
      % don't trust the calculation anyway, in particular for fractional updates...
      
      fprintf('%d: %g (stepsize: %g)\n',teller,logp,stepsize);
   
      % check whether marginal loglikelihood is going up and down, if so, lower stepsize
      
      change = logp-logpold;
      if change*oldchange < 0,   % possibly cycling
         stepsize = stepsize/2;
      end
      oldchange = change;
   end

end

if teller == opt.niter
	EP_ERR = 2;
	warning('Maximum number of iterations exceeded!\n');
end
logphist = logphist(1:teller);

comptime = toc;
fprintf('EP finished in %s seconds with EP_ERR=%d \n',num2str(comptime),EP_ERR);

%%% END MAIN


%%%%%%%%%
%
% compute the cavity approximations that result when subtracting a fraction of the term approximations 

function [myGauss,ok] = project_all(Gauss,terms,fraction)

if nargin < 3,
   fraction = 1;
end

% take out and project in moment form

% (2) cross terms between regression parameters and auxiliary parameters

[myGauss.diagC,myGauss.m] = ...
   rank_one_update(Gauss.diagC,Gauss.m,-fraction*terms.diagK,-fraction*terms.h);

 myGauss.auxC = rank_one_update(Gauss.auxC,[],-fraction*terms.auxK);

% check whether all precision matrices are strictly positive definite

if nargout > 1,
   ok =  (all(myGauss.diagC > 0) & all(myGauss.auxC > 0));
end


%%%%%%%%%
%
% compute the new term approximation by applying an adf update on all the cavity approximations

function [fullterms,crosslogz] = ...
   adfupdate_all(myGauss,terms,xhermite,whermite,xlaguerre,wlaguerre,fraction,temperature)

if nargin < 8,
   temperature = 1;
end
if nargin < 7,
   fraction = 1;
end

%% (2) cross terms between regression parameters and auxiliary variables

oldm      = myGauss.m;
oldC      = myGauss.diagC;
oldlambda = myGauss.auxC;
nfeatures = length(oldm);
nlaguerre = length(wlaguerre);

% this part heavily relies on the accompanying note
% basic idea:
% - the cavity approximation on U is an exponential distribution
% - we have analytical formulas for the moments of x conditioned upon U
% - marginal moments can then be computed through numerical integration with Gauss-Laguerre

% translate and scale the sample points to get the correct mean

U   = 2*oldlambda*xlaguerre';     % nfeatures x nlaguerre

mm = repmat(oldm,1,nlaguerre);
CC = repmat(oldC,1,nlaguerre);

% compute the partition function (integral over x) given U and turn this into weights required for computing the marginal moments

g           = (1-fraction)*log(U)/2 - fraction*mm.^2./(U + fraction*CC)/2 - log(U+ fraction*CC)/2 - (1+2*fraction)*log(2*pi)/2;
g           = g + log(repmat(wlaguerre',nfeatures,1));
maxg        = max(g,[],2);
g           = g-repmat(maxg,1,nlaguerre);
expg        = exp(g);
denominator = sum(expg,2);
neww        = expg./repmat(denominator,1,nlaguerre);

% compute the marginal moments through numerical integration

ExgU = mm.*U./(U + fraction*CC);
Ex   = sum(ExgU.*neww,2);

ExxgU = ExgU.^2 + CC.*U./(U+fraction*CC);
Exx   = sum(ExxgU.*neww,2);
EU    = sum(U.*neww,2);

newm      = Ex;
newC      = Exx-Ex.^2;
newlambda = EU/2;

% derive the term approximation from the change in mean and variance

[fullterms.diagK,fullterms.h,logzextra1] = ...
   compute_termproxy(newC,newm,oldC,oldm,fraction);

% same for the auxiliary variables, where we note that the mean will always be zero

[fullterms.auxK,dummy,logzextra2] = ...
   compute_termproxy(newlambda,zeros(nfeatures,1),oldlambda,zeros(nfeatures,1),fraction);

%crosslogz = maxg + log(denominator) + logzextra1 + (logzextra2 + logzextra2); 
crosslogz  = maxg + log(denominator);



%%%%%%%%%%
%
% compute the moments corresponding to the canonical parameters

function [Gauss,logp] = canonical_to_moments(Gauss)

[nsamples,nfeatures] = size(Gauss.A);

%% (1) regression parameters

if nsamples > nfeatures,   % in the non-degenerate case, this direct route is more stable and faster
  
   scaledA     = Gauss.A.*(repmat(Gauss.hatK,1,nfeatures));
   K           = Gauss.A'*scaledA + diag(Gauss.diagK);
   [C,logdet1] = invert_chol(K);
   Gauss.m     = C*Gauss.h;
   Gauss.diagC = diag(C);
 
else
   % this part heavily relies on the appendix of the accompanying note
   % basic idea:
   % - the precision matrix K is of the form A' hatK A + diagK, where both hatK and diagK are diagonal matrices
   % - apply Woodbury's formula to replace inverting an (nfeat x nfeat) matrix by an (nsample x nsample) alternative
   % - projections of the covariance matrix and the mean onto the feature matrix then follow immediately
   
   scaledA = Gauss.A./(repmat(Gauss.diagK',nsamples,1));
   W       = Gauss.A*scaledA';
   W       = (W + W')/2;    % make symmetric
   
   [Q,logdet1] = invert_chol(diag(1./Gauss.hatK) + W);
   logdet1     = logdet1 + sum(log(Gauss.diagK)) + sum(log(Gauss.hatK));       

   Gauss.diagC = 1./Gauss.diagK;
   for i=1:nfeatures,
      Gauss.diagC(i) = Gauss.diagC(i) - scaledA(:,i)'*Q*scaledA(:,i);
   end
   Gauss.m = Gauss.h./Gauss.diagK - scaledA'*(Q*(scaledA*Gauss.h));
end


Gauss.hatn = Gauss.A*Gauss.m;

qterm      = sum(Gauss.m .* Gauss.diagK .* Gauss.m);              % = m' * diagK * m
qterm      = qterm + sum(Gauss.hatn .* Gauss.hatK .* Gauss.hatn); % = m'*A'*K_hat*A*m
logp1      = 0.5*(qterm - logdet1);


%% (2) auxiliary variables; i.e., wrt scale mixture representation of
% Laplace prior

% this is (by far) the most expensive step when nsamples << nfeatures
% and the precision matrix of the auxiliary variables is non-diagonal

[auxC,logdet2] = invert_chol(Gauss.auxK); % only need diagonal terms
Gauss.auxC     = full(diag(auxC));        % turn into full vector
logp2          = 0.5*( - logdet2);
logp           = logp1 + 2*logp2;



%%%%%%%%%%
%
% take out the old term proxies and add the new termproxies and check whether the resulting Gaussian is still normalizable

function [newGauss,newterms,ok] = try_update(Gauss,fullterms,terms,stepsize)

if nargin < 4,
   stepsize = 1;
end

% take out the old term proxies

newGauss = update_Gauss(Gauss,terms,-1);

% compute the new term proxies as a weighted combi of the old ones and the "full" (stepsize 1) term proxies

newterms = combine_terms(fullterms,terms,stepsize);

% add the new term proxies

newGauss = update_Gauss(newGauss,newterms,1);

[L,check] = chol(newGauss.auxK,'lower');   % check whether full covariance matrix is ok
                % note that this is bit inefficient, since we redo the Cholesky later when everything is fine

ok = (check == 0 & all(newGauss.hatK > 0) & all(newGauss.diagK > 0));  % perhaps a bit too strong???


%%%%%%%%%%%%
%
% compute the moment form of the current Gauss and all cavity approximations and check whether they are fine

function [Gauss,myGauss,logdet,ok] = try_project(Gauss,terms,fraction)

if nargin < 3,
   fraction = 1;
end

[Gauss,logdet] = canonical_to_moments(Gauss);
[myGauss,ok]   = project_all(Gauss,terms,fraction);

%%%%%%%%%%
%
% if we use a temperature < 1, to get closer to the MAP solution, we have to change the prior and initial term proxies accordingly

function [Gauss,terms] = correct_for_temperature(Gauss,terms,temperature)

% note: choose temperature small to implement MAP-like behavior

Gauss.hatK = Gauss.hatK/temperature;
Gauss.h = Gauss.h/temperature;
Gauss.auxK = Gauss.auxK/temperature;
Gauss.diagK = Gauss.diagK/temperature;


% terms.hatK = terms.hatK/temperature;
% terms.hath = terms.hath/temperature;
terms.diagK = terms.diagK/temperature;
terms.auxK = terms.auxK/temperature;
terms.h = terms.h/temperature;

%%%%%%%%%%
%
% invert a positive definite matrix using Cholesky factorization

function [invA,logdet] = invert_chol(A)

[L,CHOL_ERR,S] = chol(sparse(A),'lower');   % now A = S*(L*L')*S' and (L*L') = S'*A*S
if CHOL_ERR
	error('Matrix not p.d.!');
end

invA = fastinvc(L);
invA = S*invA*S';

if nargout > 1,
   logdet = 2*sum(log(full(diag(L))));
end


%%%%%%%%%%
%
% compute the term proxy when [oldC,oldm] changes to [newC,newm]

function [K,h,logz] = compute_termproxy(newC,newm,oldC,oldm,fraction)

if nargin < 5,
   fraction = 1;
end

K = (1./newC - 1./oldC)/fraction;
h = (newm./newC - oldm./oldC)/fraction;

logz =  - log(newC./oldC)/2 + oldm.^2./oldC/2 - newm.^2./newC/2;


%%%%%%%%%%
%
% Sherman-Morrison formula to compute the change from [oldC,oldm] to [newC,newm] when we add [K,h] to the corresponding canonical parameters

function [newC,newm] = rank_one_update(oldC,oldm,K,h)

dummy = K.*oldC;
oneminusdelta = 1./(1+dummy);

newC = oneminusdelta.*oldC;

if nargout > 1,
   newm = oneminusdelta.*(oldm + h.*oldC);
end


%%%%%%%%%%%
%
% general procedure for a weighted combi of the fields of two structures

function terms = combine_terms(terms1,terms2,stepsize)

names1 = fieldnames(terms1);
names2 = fieldnames(terms2);
names = intersect(names1,names2);

terms = struct;
for i=1:length(names)   
  terms.(names{i}) = stepsize*terms1.(names{i}) + (1-stepsize)*terms2.(names{i});
end


%%%%%%%%%%%
%
% updates the Gaussian representation with new term proxies

function Gauss = update_Gauss(Gauss,terms,const)

if nargin < 3,
   const = 1;
end

Gauss.h     = Gauss.h  + const*terms.h;
Gauss.diagK = Gauss.diagK + const*terms.diagK;

% get diagonal elements
diagidx = 1:(size(Gauss.auxK,1)+1):numel(Gauss.auxK);
Gauss.auxK(diagidx) = Gauss.auxK(diagidx) + const*terms.auxK';

