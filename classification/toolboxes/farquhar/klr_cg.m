function [wb,f,J]=klr_cg(K,Y,C,varargin);
% Regularised Kernel Logistic Regression Classifier
%
% [alphab,J]=klr_cg(K,Y,lambda,varargin)
% Regularised Kernel Logistic Regression Classifier using a pre-conditioned
% conjugate gradient solver so extends to large input kernels.
%
% J = C(1) w' K w + sum_i log( (1 + exp( - y_i ( w'*K_i + b ) ) )^-1 ) 
%
% Inputs:
%  K       - NxN kernel matrix
%  Y       - Nx1 matrix of +1,-1 labels
%  C       - the regularisation parameter
%            good default is = .1*(mean(diag(K))-mean(K(:))))
%
% Outputs:
%  alphab  - (N+1)x1 matrix of the kernel weights and the bias [alpha;b]
%  p       - [Nx1] vector of conditional probabilities, Pr(y|x)
%  J       - the final objective value
%
% Options:
%  alphab  - initial guess at the kernel parameters, [alpha;b]
%  nobias  - flag we don't want the bias computed (x2 faster!)
%  maxIter - max number of CG steps to do
%  tol     - absolute error tolerance
%  tol0    - relative error tolerance, w.r.t. initial gradient
%  verb    - verbosity
%  step    - initial step size guess
%  wght    - point weights [Nx1] vector of label accuracy probabilities
%            [2x1] for per class weightings
%

% 
%
if ( nargin < 3 ) C(1)=0; end;
opts=struct('alphab',[],'nobias',0,...
            'maxIter',inf,'maxEval',10000,'tol',1e-6,'tol0',0,'objTol',0,...
            'verb',0,'step',0,'wght',[],'X',[],'ridge',1e-9,...
            'maxStep',3,'minStep',5e-2);
[opts,varargin]=parseOpts(opts,varargin{:});
opts.ridge=opts.ridge(:);

[dim,N]=size(K);

wb=opts.alphab; if ( isempty(wb) ) wb=zeros(N+1,1); end 
if ( opts.nobias ) wb(end)=0; end;

% Ensure all inputs have a consistent precision
if(isa(K,'double') & isa(Y,'single') ) Y=double(Y); end;
if(isa(K,'single')) eps=1e-7; else eps=1e-16; end;

% N.B. this form of loss weighting has no true probabilistic interpertation!
wghtY=Y;
if ( ~isempty(opts.wght) ) % point weighting -- only needed in wghtY
   if ( numel(opts.wght)==2 ) % per class weights
      wghtY(Y<0)=-1*opts.wght(1); wghtY(Y>0)=1*opts.wght(2);
   elseif ( numel(opts.wght)==N )
      wghtY=Y.*opts.wght;
   else
      error('Weight must be 2 or N elements long');
   end
end

wK   = wb(1:end-1)'*K + opts.ridge'.*wb(1:end-1)';%include ridge just in case
g    = 1./(1+exp(-Y'.*(wK+wb(end)))); g=max(g,eps); % =Pr(x|y), stop log 0
Yerr = wghtY'.*(1-g);

% precond'd gradient K^-1 (lambda*wK-K((1-g).Y)) = lambda w - (1-g).Y
MdJ   = [(2*C(1)*wb(1:end-1) - Yerr'); ...
         -sum(Yerr)];
dJ    = [K*MdJ(1:end-1)+opts.ridge*MdJ(1:end-1); ...
         MdJ(end)];
if ( opts.nobias ) MdJ(end)=0; dJ(end)=0; end;
% MdJ   = [(C(1)*wb(1:end-1) - Yerr')./diag(K); -sum(Yerr) ];
% dJ    = [K*(MdJ(1:end-1).*diag(K)); MdJ(end)];
Mr   =-MdJ;
d    = Mr;
dtdJ =-d'*dJ;
r2   = dtdJ;
r02  = r2;

Ed   = -log(g)*(Y.*wghtY); % -ln P(D|w,b,fp)
Ew   = wK*wb(1:end-1);     % -ln P(w,b|R);
J    = Ed + C(1)*Ew;       % J=neg log posterior

step=opts.step;
if( step<=0 ) step=abs(J/dtdJ); end  % init step assuming opt is at 0
step=abs(step); tstep=step;

neval=1;
if(opts.verb>0)   % debug code      
   fprintf(['%3d) %3d x=[%8g,%8g,.] J=%5g |dJ|=%8g\n'],0,neval,wb(1),wb(2),J,r2);
end

% pre-cond non-lin CG iteration
for i=1:opts.maxIter;

   oJ= J; oMr  = Mr; or2=r2; % record info about prev result we need

   %---------------------------------------------------------------------
   % Secant method for the root search.
   if ( opts.verb > 1 )
      fprintf('.%d %g=%g @ %g\n',0,0,dtdJ,J);
      if ( opts.verb>2 ) 
         hold off;plot(0,dtdJ,'r*');hold on;text(0,double(dtdJ),num2str(0)); 
         grid on;
      end
   end;
   step=max(tstep,1e-6/dtdJ); % prev step size is first guess!
   odtdJ=dtdJ; % one step before is same as current
   wb = wb + step*d;
   Kd = [K*d(1:end-1)+opts.ridge'*d(1:end-1);d(end)];%cache, so don't comp dJ
   for j=1:50;
      neval=neval+1;
      oodtdJ=odtdJ; odtdJ=dtdJ; % prev and 1 before grad values
      
      % Eval the gradient at this point.  N.B. only gradient needed for secant
      wK   = wb(1:end-1)'*K + opts.ridge'.*wb(1:end-1)'; % include ridge
      g     = 1./(1+exp(-Y'.*(wK+wb(end)))); g=max(g,eps); % stop log 0
      Yerr  = wghtY'.*(1-g);
      MdJ   = [2*C(1)*wb(1:end-1) - Yerr';...
               -sum(Yerr)];
      if ( opts.nobias ) MdJ(end)=0; end
      dtdJ   =-Kd'*MdJ;  % gradient along the line      
      % convergence test
      if( abs(dtdJ)<opts.tol | abs(dtdJ*step)<eps ) break; end; 
      
      if ( opts.verb > 1 )
         Ed   = -log(g)*(Y.*wghtY);         % P(D|w,b,fp)
         Ew   = wK*wb(1:end-1);             % P(w,b|R);
         J    = Ed + C(1)*Ew;               % J=neg log posterior         
         fprintf('.%d %g=%g @ %g\n',j,tstep,dtdJ,J); 
         if ( opts.verb > 2 ) 
            plot(tstep,dtdJ,'*'); text(double(tstep),double(dtdJ),num2str(j));
         end
      end;
      
      % convergence test, and numerical res test
      if ( abs(dtdJ) <= opts.tol || abs(dtdJ*step)<eps ) break; end; 

      % now compute the new step size
      % backeting check, so it always decreases
      if ( oodtdJ*odtdJ < 0 & odtdJ*dtdJ > 0 ...      % oodtdJ still brackets
           & abs(step*odtdJ/(odtdJ-dtdJ)) > abs(ostep) ) % would jump outside 
         step = ostep + step; % make as if we jumped here directly.
         odtdJ = oodtdJ;
      end
      ostep = step;
      % *RELATIVE* secant step size
      nstep = dtdJ/(odtdJ-dtdJ);
      nstep = sign(nstep)*max(opts.minStep,min(abs(nstep),opts.maxStep));
      step  = step * nstep ;
      tstep = tstep + step;            % total step size
      
      % move to the new point
      wb    = wb + step*d ;
   end
   if ( opts.verb > 1 ) fprintf('\n'); end;
   % compute the other bits needed for CG iteration
   dJ   = [K*MdJ(1:end-1)+opts.ridge.*MdJ(1:end-1); ...
           MdJ(end)];
   if ( opts.nobias ) dJ(end)=0; end;
   Mr =-MdJ;
   r2 =abs(Mr'*dJ); 

   % compute the function evaluation
   Ed   = -log(g)*(Y.*wghtY);         % P(D|w,b,fp)
   Ew   = wK*wb(1:end-1);             % P(w,b|R);
   J    = Ed + C(1)*Ew;               % J=neg log posterior
   if(opts.verb>0)   % debug code      
      fprintf(['%3d) %3d x=[%8g,%8g,.] J=%5g |dJ|=%8g\n'],...
              i,neval,wb(1),wb(2),J,r2);
   end   

   if ( r2<=opts.tol || r2<=r02*opts.tol0 || ...
        neval > opts.maxEval || abs(oJ-J) <= opts.objTol ...
        || oJ-J < -1e-2) % term tests
      break;
   end; 

   delta = max((Mr-oMr)'*(-dJ)/or2,0); % Polak-Ribier
   %delta = max(r2/or2,0); % Fletcher-Reeves
   d     = Mr+delta*d;     % conj grad direction
   dtdJ   =-d'*dJ;         % new search dir grad.
   if( dtdJ <= 0 )         % non-descent dir switch to steepest
      if ( opts.verb >= 1 ) fprintf('non-descent dir\n'); end;      
      d=Mr; dtdJ=-d'*dJ; 
   end; 
   
end;
if ( opts.verb >= 0 ) 
   fprintf(['%3d) %3d x=[%8g,%8g,.] J=%5g |dJ|=%8g\n'],...
           i,neval,wb(1),wb(2),J,r2);
end

% compute final decision values.
f = wb(1:end-1)'*K + wb(end);
p = 1./(1+exp(-f)); % Pr(y==1|x,w,b)
return;

%-----------------------------------------------------------------------
function [opts,varargin]=parseOpts(opts,varargin)
% refined and simplified option parser with structure flatten
i=1;
while i<=numel(varargin);  
   if ( iscell(varargin{i}) ) % flatten cells
      varargin={varargin{1:i} varargin{i}{:} varargin{i+1:end}};
   elseif ( isstruct(varargin{i}) )% flatten structures
      cellver=[fieldnames(varargin{i})'; struct2cell(varargin{i})'];
      varargin={varargin{1:i} cellver{:} varargin{i+1:end} };
   elseif( isfield(opts,varargin{i}) ) % assign fields
      opts.(varargin{i})=varargin{i+1}; i=i+1;
   else
      error('Unrecognised option');
   end
   i=i+1;
end
return;

%-----------------------------------------------------------------------------
function []=testCase()
%Make a Gaussian balls + outliers test case
[X,Y]=mkMultiClassTst([-1 0; 1 0; .2 .5],[400 400 50],[.3 .3; .3 .3; .2 .2],[],[-1 1 1]);[dim,N]=size(X);

K=X'*X; % N.B. add 1 to give implicit bias term
fInds=gennFold(Y,10,'perm',1); trnInd=any(fInds(:,1:9),2); tstInd=fInds(:,10);
trnSet=find(trnInd);

[alphab,J]=klr(K(trnInd,trnInd),Y(trnInd),1,'verb',1);
dv=K(tstInd,trnInd)*alphab(1:end-1)+alphab(end);
dv2conf(Y(tstInd),dv)

% for linear kernel
alpha=zeros(N,1);alpha(find(trnInd))=alphab(1:end-1); % equiv alpha
plotLinDecisFn(X,Y,X(:,trnInd)*alphab(1:end-1),alphab(end),alpha);

% unbalanced data
wght=[1,sum(Y>0)/sum(Y<=0)];
