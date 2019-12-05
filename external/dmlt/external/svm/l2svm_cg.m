function [wb,f,J,obj]=l2svm_cg(K,Y,C,varargin)
% [alphab,f,J]=l2svm(K,Y,C,varargin)
% Quadratic Loss Support Vector machine using a pre-conditioned conjugate
% gradient solver so extends to large input kernels.
%
% J = C(1) w' K w + sum_i max(0 , 1 - y_i ( w'*K_i + b ) ).^2
%
% Inputs:
%  K       - [NxN] kernel matrix
%  Y       - [Nx1] matrix of -1/0/+1 labels, (0 label pts are implicitly ignored)
%  C       - the regularisation parameter, roughly max allowed length of the weight vector
%            good default is: .1*var(data) = .1*(mean(diag(K))-mean(K(:))))
%
% Outputs:
%  alphab  - [(N+1)x1] matrix of the kernel weights and the bias [alpha;b]
%  f       - [Nx1] vector of decision values
%  J       - the final objective value
%  obj     - [J Ed Ew]
%  p       - [Nx1] vector of conditional probabilities, Pr(y|x)
%
% Options:
%  alphab  - [(N+1)x1] initial guess at the kernel parameters, [alpha;b] ([])
%  ridge   - [float] ridge to add to the kernel to improve convergence.
%             ridge<0 -- absolute ridge value
%             ridge>0 -- size relative to the mean kernel eigenvalue
%  maxEval - [int] max number for function evaluations                    (N*5)
%  maxIter - [int] max number of CG steps to do                           (inf)
%  maxLineSrch - [int] max number of line search iterations to perform    (50)
%  objTol0 - [float] relative objective gradient tolerance                (1e-5)
%  objTol  - [float] absolute objective gradient tolerance                (0)
%  tol0    - [float] relative gradient tolerance, w.r.t. initial value    (0)
%  lstol0  - [float] line-search relative gradient tolerance, w.r.t. initial value   (1e-2)
%  tol     - [float] absolute gradient tolerance                          (0)
%  verb    - [int] verbosity                                              (0)
%  step    - initial step size guess                                      (1)
%  wght    - point weights [Nx1] vector of label accuracy probabilities   ([])
%            [2x1] for per class weightings
%            [1x1] relative weight of the positive class
%  nobias  - [bool] flag we don't want the bias computed                  (false)
%  issqrtmK - [bool] flag to indicate whether the K represents the possibly
%                    low dimensional sqrtm of K, for computational speed
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)

% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied
if nargin < 3
  C(1)=0;
end
opts = struct('alphab',  [],  ...
              'nobias',  0,   ...
              'dim',     [],  ...
              'maxIter', inf, ...
              'maxEval', [],  ...
              'tol',     0,   ...
              'tol0',    0,   ... 
              'lstol0',  1e-4,...
              'objTol',  0,   ...
              'objTol0', 1e-5,...
              'verb',    0,   ...
              'step',    0,   ...
              'wght',    [],  ...
              'X',       [],  ...
              'ridge',   0,   ...
              'maxLineSrch', 50, ...
              'maxStep', 3,   ...
              'minStep', 5e-2,...
              'weightDecay', 0, ...
              'marate',  .95, ...
              'bPC',     [],  ...
              'incThresh', .75, ...
              'issqrtmK', false);
opts       = parseOpts(opts,varargin{:});
opts.ridge = opts.ridge(:);
if isempty(opts.maxEval)
  opts.maxEval = 5*sum(Y(:)~=0); 
end

% Ensure all inputs have a consistent precision
if isa(K,'double') && isa(Y,'single')
  Y = double(Y); 
end
if isa(K,'single')
  eps=1e-7; 
else
  eps=1e-16;
end
opts.tol = max(opts.tol,eps); % gradient magnitude tolerence

N  = size(K,1);
Y  = Y(:); % ensure Y is col vector

% check for degenerate inputs
if all(Y>=0) || all(Y<=0)
  warning('Degnerate inputs, 1 class problem');
end

if opts.ridge>0 % make the ridge relative to the max eigen-value
  opts.ridge = opts.ridge*median(abs(diag(K)));
  ridge = opts.ridge;
else % negative value means absolute ridge
  ridge = abs(opts.ridge);
end

% generate an initial seed solution if needed
wb = opts.alphab;   % N.B. set the initial solution
if isempty(wb)
  wb = zeros(N+1,1,class(K));
  
  % prototype classifier seed
  wb(Y>0) =  .5./sum(Y>0); 
  wb(Y<0) = -.5./sum(Y<0); % vector between pos/neg class centers
  
  if opts.issqrtmK
    wK = (wb(1:end-1)'*K)*K';
  else
    wK = wb(1:end-1)'*K;
  end
  wb(end) = -(wK(Y<0)*abs(wb(Y<0))*sum(Y>0)  + wK(Y>0)*abs(wb(Y>0))*sum(Y<0))./sum(Y~=0); % weighted mean
  
  % find least squares optimal scaling and bias
  sb = [C(1)*wK*wb(1:end-1)+wK*wK' sum(wK); sum(wK) sum(Y~=0)]\[wK*Y; sum(Y)];
  
  wb(1:end-1) = wb(1:end-1)*sb(1); 
  wb(end)     = sb(2);
end
if opts.nobias
  wb(end) = 0;
end;

% check if it's more efficient to sub-set the kernel, because of lots of ignored points
oK = K;
oY = Y;
incIdx = (Y(:)~=0);
if sum(incIdx)/numel(Y) < opts.incThresh % if enough ignored to be worth it
  if sum(incIdx)==0
    error('Empty training set!');
  end
  K  = K(incIdx,incIdx);
  Y  = Y(incIdx);
  wb = wb([incIdx; true]);
end

wght  = 1; 
wghtY = Y;
if ~isempty(opts.wght) % point weighting -- only needed in wghtY
  if numel(opts.wght)==1 % weight ratio between classes
    wght      = zeros(size(Y));
    wght(Y<0) =  1./sum(Y<0); 
    wght(Y>0) = (1./sum(Y>0))*opts.wght;
    wght      = wght*sum(abs(Y))./sum(abs(wght)); % ensure total weighting is unchanged
  elseif numel(opts.wght)==2 % per class weights
    wght      = zeros(size(Y));
    wght(Y<0) = 1*opts.wght(1); 
    wght(Y>0) = 1*opts.wght(2);
  elseif numel(opts.wght)==N
  else
    error('Weight must be 2 or N elements long');
  end
  wghtY = wght.*Y;
else
  wght = 1;
end

% Normalise the kernel to prevent rounding issues causing convergence problems
% = average kernel eigen-value + regularisation const = ave row norm
if opts.issqrtmK
  diagK = sum(K.^2,2);
else
  diagK = K(1:size(K,1)+1:end);
end
if sum(incIdx)<size(K,1) 
  diagK = diagK(incIdx); 
end
muEig = median(diagK); % approx hessian scaling, for numerical precision
% adjust alpha and regul-constant to leave solution unchanged
wb(1:end-1) = wb(1:end-1)*muEig;
C(1)        = C(1)./muEig;

% set the bias (i.e. b) pre-conditioner
bPC = opts.bPC;
if ( isempty(bPC) ) % bias pre-condn with the diagonal of the hessian
  bPC  = sqrt(abs(muEig + 2*C(1))./muEig);   % N.B. use sqrt for safety?
  bPC  = 1./bPC;
end

% include ridge and re-scale by muEig for numerical stability
if opts.issqrtmK
  wK = ( (wb(1:end-1)'*K)*K' + ridge'.*wb(1:end-1)')./muEig; % include ridge
else
  wK = (wb(1:end-1)'*K + ridge'.*wb(1:end-1)')./muEig; % include ridge
end
err  = 1-Y.*(wK'+wb(end));
svs  = err>0 & Y~=0;
Yerr = wghtY.*err;  % weighted error
% pre-conditioned gradient,
% K^-1*dJdw = K^-1(2 C(1)Kw - 2 K I_sv(Y-f)) = 2*(C(1)w - Isv (Y-f) )
MdJ  = [(2*C(1)*wb(1:end-1) - 2*(Yerr.*svs)); -2*sum(Yerr(svs))./bPC];
if opts.issqrtmK
  dJ = [(K*(K'*MdJ(1:end-1))+ridge*MdJ(1:end-1))./muEig; -2*sum(Yerr(svs))];
else
  dJ = [(K*MdJ(1:end-1)+ridge*MdJ(1:end-1))./muEig; -2*sum(Yerr(svs))];
end
if opts.nobias 
  MdJ(end) = 0; 
  dJ(end)  = 0; 
end
Mr   =-MdJ;
d    = Mr;
dtdJ =-(d'*dJ);
r2   = dtdJ;
r02  = r2;

Ed   = (wght.*err(svs))'*err(svs);
Ew   = wK*wb(1:end-1);
J    = Ed + C(1)*Ew; % SVM objective

% Set the initial line-search step size
step = opts.step;
if step<=0
  step = min(sqrt(abs(J/max(dtdJ,eps))),1); 
end %init step assuming opt is at 0
step  = abs(step); 
tstep = step;

neval = 1;
lend  = '\r';
if opts.verb>0    % debug code
  if opts.verb>1
    lend='\n';
  else
    fprintf('\n');
  end;
  fprintf('%3d) %3d x=[%5f,%5f,.] J=%5f (%5f+%5f) |dJ|=%8g\n',0,neval,wb(1),wb(2),J,Ew./muEig,Ed,r2);
end

% pre-cond non-lin CG iteration
J0   = J; 
madJ = abs(J); % init-grad est is init val
wb0  = wb; 
Kd   = zeros(size(wb),class(wb)); 
dJ   = zeros(size(wb),class(wb));
for iter = 1:min(opts.maxIter,2e6)  % stop some matlab versions complaining about index too big
  
  oJ  = J; 
  oMr = Mr; 
  or2 = r2; 
  owb = wb; % record info about prev result we need
  
  %---------------------------------------------------------------------
  % Secant method for the root search.
  if opts.verb > 2
    fprintf('.%d %g=%g @ %g (%g+%g)\n',0,0,dtdJ,J,Ed,Ew./muEig);
    if ( opts.verb>3 )
      hold off;plot(0,dtdJ,'r*');hold on;text(0,double(dtdJ),num2str(0));
      grid on;
    end
  end
  ostep = inf;
  step  = tstep; % max(tstep,abs(1e-6/dtdJ)); % prev step size is first guess!
  odtdJ = dtdJ;  % one step before is same as current
  
  wK0 = wK;
  if opts.issqrtmK
    dK  = ( (d(1:end-1)'*K)*K'+d(1:end-1)'.*ridge')./muEig;
  else
    dK  = (d(1:end-1)'*K+d(1:end-1)'.*ridge')./muEig; % N.B. v'*M is 50% faster than M*v'!!!
  end  
  db  = d(end);
  dKw = dK*wb(1:end-1); 
  dKd = dK*d(1:end-1);
  
  dtdJ0 = abs(dtdJ); % initial gradient, for Wolfe 2 convergence test
  for j=1:opts.maxLineSrch
    neval  = neval+1;
    oodtdJ = odtdJ; 
    odtdJ  = dtdJ; % prev and 1 before grad values
    
    % Eval the gradient at this point.  N.B. only gradient needed for secant
    wK   = wK0 + tstep*dK;
    err  = 1-Y.*(wK'+wb(end)+tstep*d(end)); 
    svs  = err>0 & Y~=0;
    Yerr = wghtY.*err;
    dtdJ = -(2*C(1)*(dKw+tstep*dKd) - 2*dK*(Yerr.*svs) + -2*db*sum(Yerr(svs))); % gradient along the line @ new position
    
    if opts.verb > 2
      Ed     = (wght.*err(svs))*err(svs)';
      Ew     = wK*wb(1:end-1);
      J      = Ed + C(1)*Ew; % SVM objective
      fprintf('.%d %g=%g @ %g (%g+%g)\n',j,tstep,dtdJ,J,Ed,Ew./muEig);
      if opts.verb > 3
        plot(tstep,dtdJ,'*'); text(double(tstep),double(dtdJ),num2str(j));
      end
    end
    
    % convergence test, and numerical res test
    if iter>1 || j>3 % Ensure we do decent line search for 1st step size!
      if abs(dtdJ) < opts.lstol0*abs(dtdJ0) || ... % Wolfe 2, gradient enough smaller
          abs(dtdJ*step) <= opts.tol               % numerical resolution
        break;
      end
    end
    
    % now compute the new step size
    % backeting check, so it always decreases
    if oodtdJ*odtdJ < 0 && odtdJ*dtdJ > 0 ...                 % oodtdJ still brackets
        && abs(step*dtdJ) > abs(odtdJ-dtdJ)*(abs(ostep+step)) % would jump outside
      step = ostep + step; % make as if we jumped here directly.
      % but prev points gradient, this is necessary stop very steep orginal gradient preventing decent step sizes
      odtdJ = -sign(odtdJ)*sqrt(abs(odtdJ))*sqrt(abs(oodtdJ)); % geometric mean
    end
    ostep = step;
    % *RELATIVE* secant step size
    ddtdJ = odtdJ-dtdJ;
    if ddtdJ~=0
      nstep = dtdJ/ddtdJ;
    end % secant step size, guard div by 0
    nstep = sign(nstep)*max(opts.minStep,min(abs(nstep),opts.maxStep)); % bound growth/min-step size
    step  = step * nstep ;           % absolute step
    tstep = tstep + step;            % total step size
    
    % % move to the new point
    % wb    = wb + step*d ;
  end
  if opts.verb > 2
    fprintf('\n'); 
  end
  % update the solution with this step
  wb  = wb + tstep*d;
  
  % compute the other bits needed for CG iteration
  MdJ  = [(2*C(1)*wb(1:end-1) - 2*(Yerr.*svs)); -2*sum(Yerr(svs))./bPC];
  if opts.issqrtmK
    dJ(1:end-1) = ((MdJ(1:end-1)'*K)*K'+MdJ(1:end-1)'.*ridge)./muEig;
  else
    dJ(1:end-1) = (MdJ(1:end-1)'*K+MdJ(1:end-1)'.*ridge)./muEig;
  end
  dJ(end)     = bPC*MdJ(end);
  if opts.nobias
    dJ(end)=0;
  end
  Mr = -MdJ;
  r2 = abs(Mr'*dJ);
  
  % compute the function evaluation
  Ed  = (wght.*err(svs))'*err(svs);
  Ew  = wK*wb(1:end-1);
  J   = Ed + C(1)*Ew; % SVM objective
  if(opts.verb>0)   % debug code
    fprintf(['%3d) %3d x=[%8f,%8f,.] J=%5f (%5f+%5f) |dJ|=%8g' lend],...
      iter,neval,wb(1),wb(2),J,Ew./muEig,Ed,r2);
  end
  
  if J > oJ*(1+1e-3) || isnan(J) % check for stuckness
    if opts.verb>=1 
      warning('Line-search Non-reduction - aborted'); 
    end
    J  = oJ; 
    wb = owb; 
    break
  end
  
  %------------------------------------------------
  % convergence test
  if iter==1
    madJ = abs(oJ-J);
    dJ0  = max(abs(madJ),eps); 
    r02  = r2;
  elseif iter<5
    dJ0  = max(dJ0,abs(oJ-J));
    r02  = max(r02,r2); % conv if smaller than best single step
  end
  madJ = madJ*(1-opts.marate)+abs(oJ-J)*(opts.marate); %move-ave objective grad est
  if r2<=opts.tol || ... % small gradient + numerical precision
     r2< r02*opts.tol0 || ... % Wolfe condn 2, gradient enough smaller
     neval > opts.maxEval || ... % abs(odtdJ-dtdJ) < eps || ... % numerical resolution
     madJ <= opts.objTol || madJ < opts.objTol0*dJ0 % objective function change
    break
  end
  
  %------------------------------------------------
  % conjugate direction selection
  delta = max((Mr-oMr)'*(-dJ)/or2,0); % Polak-Ribier
  %delta = max(r2/or2,0); % Fletcher-Reeves
  d     = Mr+delta*d;     % conj grad direction
  dtdJ  = -d'*dJ;         % new search dir grad.
  if dtdJ <= 0         % non-descent dir switch to steepest
    if opts.verb >= 2
      fprintf('non-descent dir\n'); 
    end
    d    = Mr;
    dtdJ = -d'*dJ;
  end
  
  if opts.weightDecay > 0 % decay term to 0 non-svs faster
    pts = ~svs' & wb(1:end-1).*d(1:end-1) < 0;
    %wb(pts)=wb(pts)*opts.weightDecay;
    d(pts) = d(pts)*opts.weightDecay;
  end
end

if opts.verb >= 0
  fprintf(['%3d) %3d x=[%8f,%8f,.] J=%5f (%5f+%5f) |dJ|=%8g\n'],...
    iter,neval,wb(1),wb(2),J,Ew./muEig,Ed,r2);
end

if J > J0*(1+1e-4) || isnan(J)
  if opts.verb>=0
    warning('Non-reduction'); 
  end
  wb = wb0;
end

% fix the stabilising K normalisation
wb(1:end-1) = wb(1:end-1)./muEig;

% compute final decision values.
if numel(Y)~=numel(incIdx) % map back to the full kernel space, if needed
  nwb = zeros(size(oK,1)+1,1); 
  nwb(incIdx) = wb(1:end-1); 
  nwb(end)    = wb(end); 
  wb  = nwb;
  K   = oK;
  Y   = oY;
end

if opts.issqrtmK
  f = (wb(1:end-1)'*K)*K' + wb(end); f = reshape(f,size(Y));
else
  f = wb(1:end-1)'*K + wb(end); f = reshape(f,size(Y));
end
obj = [J Ew./muEig Ed];

%-----------------------------------------------------------------------
function [opts,varargin] = parseOpts(opts,varargin)
% refined and simplified option parser with structure flatten
i = 1;
while i<=numel(varargin)
  if iscell(varargin{i}) % flatten cells
    varargin = {varargin{1:i} varargin{i}{:} varargin{i+1:end}};
  elseif isstruct(varargin{i})% flatten structures
    cellver  = [fieldnames(varargin{i})'; struct2cell(varargin{i})'];
    varargin = {varargin{1:i} cellver{:} varargin{i+1:end} };
  elseif isfield(opts,varargin{i}) % assign fields
    opts.(varargin{i})=varargin{i+1}; 
    i = i+1;
  else
    error('Unrecognised option');
  end
  i = i+1;
end

%-----------------------------------------------------------------------------
function []=testCase()
%TESTCASE 1) hinge vs. logistic (unregularised)
[X,Y]=mkMultiClassTst([-1 0; 1 0; .2 .5],[400 400 50],[.3 .3; .3 .3; .2 .2],[],[-1 1 1]);


K=X'*X; % Simple linear kernel
N=size(X,2);
fInds=gennFold(Y,10,'perm',1); trnInd=any(fInds(:,1:9),2); tstInd=fInds(:,10);
[alphab,J]=l2svm_cg(K(trnInd,trnInd),Y(trnInd),1,'verb',2);
dv=K(tstInd,trnInd)*alphab(1:end-1)+alphab(end);
dv2conf(Y(tstInd),dv)

% for linear kernel
alpha=zeros(N,1);alpha(find(trnInd))=alphab(1:end-1); % equiv alpha
w=X(:,trnInd)*alphab(1:end-1); b=alphab(end); plotLinDecisFn(X,Y,w,b,alpha);

% check the primal dual gap!
w=alphab(1:end-1); b=alphab(end);
Jd = C(1)*(2*sum(w.*Y(trnInd)) - C(1)*w'*w - w'*K(trnInd,trnInd)*w);

% imbalance test
[X,Y]=mkMultiClassTst([-1 0; 1 0],[400 400],[.5 .5; .5 .5],[],[-1 1]);

[X,Y]=mkMultiClassTst([-1 0; 1 0],[400 40],[.5 .5; .5 .5],[],[-1 1]);
[alphab,J]=l2svm_cg(K(trnInd,trnInd),Y(trnInd),1,'verb',2,'wght',[.1 1]);

