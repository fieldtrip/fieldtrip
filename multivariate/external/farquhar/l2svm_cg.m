function [wb,f,J,obj]=l2svm_cg(K,Y,C,varargin);
% [alphab,f,J]=l2svm(K,Y,C,varargin)
% Quadratic Loss Support Vector machine using a pre-conditioned conjugate
% gradient solver so extends to large input kernels.
%
% J = C(1) w' K w + sum_i max(0 , 1 - y_i ( w'*K_i + b ) ).^2 
% 
% Inputs:
%  K       - NxN kernel matrix
%  Y       - Nx1 matrix of +1,-1 labels (0 label pts are implicitly ignored)
%  C       - the regularisation parameter
%            good default is = .1*(mean(diag(K))-mean(K(:))))
%
% Outputs:
%  alphab  - (N+1)x1 matrix of the kernel weights and the bias b
%  f       - [Nx1] vector of decision values
%  J       - the final objective value
%
% Options:
%  alphab  - initial guess at the kernel parameters and bias
%  maxIter - max number of CG steps to do
%  tol     - absolute error tolerance
%  tol0    - relative error tolerance, w.r.t. initial gradient
%  verb    - verbosity
%  step    - initial step size guess
%  wght    - point weights [Nx1] vector of label accuracy probabilities
%            [2x1] for per class weightings
%
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)

% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied
if ( nargin < 3 ) C(1)=0; end;
opts=struct('alphab',[],'nobias',0,'dim',[],...
            'maxIter',inf,'maxEval',[],'tol',0,'tol0',0,'lstol0',1e-2,'objTol',1e-6,'objTol0',1e-4,...
            'verb',0,'step',0,'wght',[],'X',[],'ridge',0,'maxLineSrch',50,...
            'maxStep',3,'minStep',5e-2,'weightDecay',0,'marate',.95,'bPC',[],'incThresh',.75,'optBias',1,'maxTr',1000);
[opts,varargin]=parseOpts(opts,varargin{:});
opts.ridge=opts.ridge(:);
if ( isempty(opts.maxEval) ) opts.maxEval=5*sum(Y(:)~=0); end
% Ensure all inputs have a consistent precision
if(isa(K,'double') & isa(Y,'single') ) Y=double(Y); end;
if(isa(K,'single')) eps=1e-7; else eps=1e-16; end;
opts.tol=max(opts.tol,eps); % gradient magnitude tolerence

[dim,N]=size(K); Y=Y(:); % ensure Y is col vector

wb=opts.alphab; 
if ( isempty(wb) )    
   wb=zeros(N+1,1,class(K)); 
   % prototype classifier seed
   wb(Y>0)=.5./sum(Y>0); wb(Y<0)=-.5./sum(Y<0); % vector between pos/neg class centers
   wK = wb(1:end-1)'*K;
   %wb(end)=-(wK*abs(wb(1:end-1))); % bias removes ave result
   wb(end) = -(wK(Y<0)*abs(wb(Y<0))*sum(Y>0)  + wK(Y>0)*abs(wb(Y>0))*sum(Y<0))./sum(Y~=0); % weighted mean
   wb=wb./(mean(abs(wK(:)+wb(end)))); % set to unit norm prediction value
   wb=wb*min(1,N./max(eps,C(1)));  % re-scale to near optimal size
end 
if ( opts.nobias ) wb(end)=0; end;

if ( opts.ridge>0 ) % make the ridge relative to the max eigen-value
   opts.ridge = opts.ridge*median(abs(diag(K)));
   ridge = opts.ridge;
else % negative value means absolute ridge
   ridge = abs(opts.ridge);
end

wghtY=Y;
if ( ~isempty(opts.wght) ) % point weighting -- only needed in wghtY
   if ( numel(opts.wght)==1 ) % weight ratio between classes
      wghtY(Y<0)=-1./sum(Y<0); wghtY(Y>0)=(1./sum(Y>0))*opts.wght;
      wghtY = wghtY*sum(abs(Y))./sum(abs(wghtY)); % ensure total weighting is unchanged
   elseif ( numel(opts.wght)==2 ) % per class weights
      wghtY(Y<0)=-1*opts.wght(1); wghtY(Y>0)=1*opts.wght(2);
   elseif ( numel(opts.wght)==N )
      wghtY=Y.*opts.wght;
   else
      error('Weight must be 2 or N elements long');
   end
end

% check if it's more efficient to sub-set the kernel, because of lots of ignored points
oK=K; oY=Y;
incIdx=(Y(:)~=0);
if ( sum(incIdx)./numel(Y) < opts.incThresh ) % if enough ignored to be worth it
   if ( sum(incIdx)==0 ) error('Empty training set!'); end;
   K=K(incIdx,incIdx); Y=Y(incIdx); wghtY=wghtY(incIdx); wb=wb([incIdx; true]);
end

if ( size(K,1)>=1.5*opts.maxTr && all(wb==0) ) % use sub-set to get seed solution
   idx =false(size(Y)); 
   pIdx=find(Y>0);perm=randperm(numel(pIdx));idx(pIdx(perm(1:min(end,floor(opts.maxTr./2)))))=true;
   nIdx=find(Y<0);perm=randperm(numel(nIdx));idx(nIdx(perm(1:min(end,floor(opts.maxTr./2)))))=true;
   wb0 = klr_cg(K(idx,idx),Y(idx),C,opts);
   wb(idx) = wb0(1:end-1); wb(end)=wb0(end);
end

% Normalise the kernel to prevent rounding issues causing convergence problems
% = average kernel eigen-value + regularisation const = ave row norm
diagK= K(1:size(K,1)+1:end); 
if ( sum(incIdx)>size(K,1) ) diagK=diagK(incIdx); end;
muEig=median(diagK); % approx hessian scaling, for numerical precision
% adjust alpha and regul-constant to leave solution unchanged
wb(1:end-1)=wb(1:end-1)*muEig;
C(1) = C(1)./muEig;

% set the bias (i.e. b) pre-conditioner
bPC=opts.bPC;
if ( isempty(bPC) ) % bias pre-condn with the diagonal of the hessian
   bPC  = sqrt(abs(muEig + 2*C(1))./muEig);   % N.B. use sqrt for safety?
   bPC  = 1./bPC;
   %fprintf('bPC=%g\n',bPC);
end

wK   = wb(1:end-1)'*K + ridge'.*wb(1:end-1)';%include ridge just in case
err  = 1-Y'.*(wK+wb(end)); svs=err>=0 & Y'~=0;
% pre-condinationed gradient, 
Yerr = wghtY'.*err; % wght.*Y.*(1-Y.*f)
% K^-1*dJdw = K^-1(2 C(1)Kw - 2 K I_sv(Y-f)) = 2*(C(1)w - Isv (Y-f) )
% N.B. could include additional diag(K) pre-cond here if it would help
MdJ  = [(2*C(1)*wb(1:end-1) - 2*(Yerr.*svs)')./diag(K);...
        -2*sum(Yerr(svs))/bPC];
dJ   = [K*(MdJ(1:end-1).*diag(K))+ridge.*MdJ(1:end-1).*diag(K); ...
        -2*sum(Yerr(svs))];
if ( opts.nobias ) MdJ(end)=0; dJ(end)=0; end;
Mr   =-MdJ;
d    = Mr;
dtdJ =-(d'*dJ);
r2   = dtdJ;
r02  = r2;

Ed   = (Y(svs)'.*Yerr(svs))*err(svs)';  
Ew   = wK*wb(1:end-1);
J    = Ed + C(1)*Ew; % SVM objective

% Set the initial line-search step size
step=opts.step;
if( step<=0 ) step=min(sqrt(abs(J/max(dtdJ,eps))),1); end %init step assuming opt is at 0
step=abs(step); tstep=step;

neval=1; lend='\r';
if(opts.verb>0)   % debug code      
   if ( opts.verb>1 ) lend='\n'; else fprintf('\n'); end;
   fprintf(['%3d) %3d x=[%5f,%5f,.] J=%5f (%5f+%5f) |dJ|=%8g\n'],0,neval,wb(1),wb(2),J,Ew,Ed,r2);
end

% pre-cond non-lin CG iteration
J0=J; madJ=abs(J); % init-grad est is init val
wb0=wb;
for iter=1:min(opts.maxIter,2e6);  % stop some matlab versions complaining about index too big

   oJ= J; oMr  = Mr; or2=r2; owb=wb; % record info about prev result we need

   %---------------------------------------------------------------------
   % Secant method for the root search.
   if ( opts.verb > 2 )
      fprintf('.%d %g=%g @ %g\n',0,0,dtdJ,J);
      if ( opts.verb>3 ) 
         hold off;plot(0,dtdJ,'r*');hold on;text(0,double(dtdJ),num2str(0)); 
         grid on;
      end
   end;
   ostep=inf;step=tstep;%max(tstep,abs(1e-6/dtdJ)); % prev step size is first guess!
   odtdJ=dtdJ; % one step before is same as current
   wb = wb + step*d;
   Kd = [(K*d(1:end-1)+ridge'*d(1:end-1))./muEig;bPC*d(end)];%cache, so don't comp dJ
   dtdJ0=abs(dtdJ); % initial gradient, for Wolfe 2 convergence test
   for j=1:opts.maxLineSrch;
      neval=neval+1;
      oodtdJ=odtdJ; odtdJ=dtdJ; % prev and 1 before grad values
      
      % Eval the gradient at this point.  N.B. only gradient needed for secant
      wK   = (wb(1:end-1)'*K + ridge'.*wb(1:end-1)')./muEig; % include ridge
      err  = 1-Y'.*(wK+wb(end)); svs=err>0 & Y'~=0;
      Yerr = wghtY'.*err; 
      MdJ  = [(2*C(1)*wb(1:end-1) - 2*(Yerr.*svs)'); ...
              -2*sum(Yerr(svs))./bPC];
      if ( opts.nobias ) MdJ(end)=0; end;
      dtdJ  =-Kd'*MdJ;  % gradient along the line      
      
      if ( opts.verb > 2 )
         Ed     = err(svs)*err(svs)';  Ew     = wK*wb(1:end-1);
         J      = Ed + C(1)*Ew; % SVM objective
         fprintf('.%d %g=%g @ %g \n',j,tstep,dtdJ,J); 
         if ( opts.verb > 3 )
            plot(tstep,dtdJ,'*'); text(double(tstep),double(dtdJ),num2str(j));
         end
      end;

      % convergence test, and numerical res test
      if(iter>1|j>3) % Ensure we do decent line search for 1st step size!
         if ( abs(dtdJ) < opts.lstol0*abs(dtdJ0) | ... % Wolfe 2, gradient enough smaller
              abs(dtdJ*step) <= opts.tol )              % numerical resolution
            break;
         end
      end
            
      % now compute the new step size
      % backeting check, so it always decreases
      if ( oodtdJ*odtdJ < 0 & odtdJ*dtdJ > 0 ...      % oodtdJ still brackets
           & abs(step*dtdJ) > abs(odtdJ-dtdJ)*(abs(ostep+step)) ) % would jump outside 
         step = ostep + step; % make as if we jumped here directly.
         % but prev points gradient, this is necessary stop very steep orginal gradient preventing decent step sizes
         odtdJ = -sign(odtdJ)*sqrt(abs(odtdJ))*sqrt(abs(oodtdJ)); % geometric mean
      end
      ostep = step;
      % *RELATIVE* secant step size
      ddtdJ = odtdJ-dtdJ; 
      if ( ddtdJ~=0 ) nstep = dtdJ/ddtdJ; end; % secant step size, with guard against div by 0
      nstep = sign(nstep)*max(opts.minStep,min(abs(nstep),opts.maxStep)); % bound growth/min-step size
      step  = step * nstep ;           % absolute step
      tstep = tstep + step;            % total step size
      
      % move to the new point
      wb    = wb + step*d ;
   end
   if ( opts.verb > 1 ) fprintf('\n'); end;

   
   % compute the other bits needed for CG iteration
   dJ   = [(K*MdJ(1:end-1)+ridge.*MdJ(1:end-1))./muEig; ...
           bPC*MdJ(end)];
   if ( opts.nobias ) dJ(end)=0; end;
   Mr =-MdJ;
   r2 =abs(Mr'*dJ); 
   
   % compute the function evaluation
   Ed     = (Y(svs)'.*Yerr(svs))*err(svs)';  
   Ew     = wK*wb(1:end-1);
   J      = Ed + C(1)*Ew; % SVM objective
   if(opts.verb>0)   % debug code      
      fprintf(['%3d) %3d x=[%8f,%8f,.] J=%5f (%5f+%5f) |dJ|=%8g' lend],...
              iter,neval,wb(1),wb(2),J,Ew,Ed,r2);
   end   

   if ( J > oJ*(1+1e-3) || isnan(J) ) % check for stuckness
      if ( opts.verb>=0 ) warning('Line-search Non-reduction - aborted'); end;
      J=oJ; wb=owb; break;
   end;
   
   %------------------------------------------------
   % convergence test
   if ( iter==1 )     madJ=abs(oJ-J); dJ0=max(abs(madJ),eps); r02=r2;
   elseif( iter<5 )   dJ0=max(dJ0,abs(oJ-J)); r02=max(r02,r2); % conv if smaller than best single step
   end
   madJ=madJ*(1-opts.marate)+abs(oJ-J)*(opts.marate);%move-ave objective grad est
   if ( r2<=opts.tol || ... % small gradient + numerical precision
        r2< r02*opts.tol0 || ... % Wolfe condn 2, gradient enough smaller
        neval > opts.maxEval || ... % abs(odtdJ-dtdJ) < eps || ... % numerical resolution
        madJ <= opts.objTol || madJ < opts.objTol0*dJ0 ) % objective function change
      break;
   end;    
   
   %------------------------------------------------
   % conjugate direction selection
   delta = max((Mr-oMr)'*(-dJ)/or2,0); % Polak-Ribier
   %delta = max(r2/or2,0); % Fletcher-Reeves
   d     = Mr+delta*d;     % conj grad direction
   dtdJ   =-d'*dJ;         % new search dir grad.
   if( dtdJ <= 0 )         % non-descent dir switch to steepest
      if ( opts.verb >= 2 ) fprintf('non-descent dir\n'); end;      
      d=Mr; dtdJ=-d'*dJ; 
   end; 

   if ( opts.weightDecay > 0 ) % decay term to 0 non-svs faster
      pts=~svs' & wb(1:end-1).*d(1:end-1) < 0;
      %wb(pts)=wb(pts)*opts.weightDecay; 
      d(pts) =d(pts)*opts.weightDecay;
   end
  
end;
if ( opts.verb >= 0 ) 
   fprintf(['%3d) %3d x=[%8f,%8f,.] J=%5f (%5f+%5f) |dJ|=%8g\n'],...
           iter,neval,wb(1),wb(2),J,Ew,Ed,r2);
end

if ( J > J0*(1+1e-4) || isnan(J) ) 
   if ( opts.verb>=0 ) warning('Non-reduction');  end;
   wb=wb0;
end;

% fix the stabilising K normalisation
wb(1:end-1) = wb(1:end-1)./muEig;

% compute final decision values.
if ( numel(Y)~=numel(incIdx) ) % map back to the full kernel space, if needed
   nwb=zeros(size(oK,1)+1,1); nwb(incIdx)=wb(1:end-1); nwb(end)=wb(end); wb=nwb;
   K=oK; Y=oY;
end

f = wb(1:end-1)'*K + wb(end); f = reshape(f,size(Y));
obj = [J Ew Ed];
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

