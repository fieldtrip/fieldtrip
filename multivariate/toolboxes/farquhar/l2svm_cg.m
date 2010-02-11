function [wb,f,J]=l2svm_cg(K,Y,C,varargin);
% [alphab,f,J]=l2svm(K,Y,C,varargin)
% Quadratic Loss Support Vector machine using a pre-conditioned conjugate
% gradient solver so extends to large input kernels.
%
% J = C(1) w' K w + sum_i max(0 , 1 - y_i ( w'*K_i + b ) ).^2 
% 
% Inputs:
%  K       - NxN kernel matrix
%  Y       - Nx1 matrix of +1,-1 labels
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

if ( nargin < 3 ) C=0; end;
opts=struct('alphab',[],...
            'maxIter',inf,'maxEval',inf,'tol',1e-6,'tol0',0,'objTol',0,...
            'verb',0,'step',0,'wght',[],'X',[],'ridge',1e-9,'nobias',0,...
            'maxStep',4,'minStep',1e-3,'weightDecay',0);
[opts,varargin]=parseOpts(opts,varargin{:});
opts.ridge=opts.ridge(:);

[N,dim]=size(K);

% Ensure all inputs have a consistent precision
if(isa(K,'double') & isa(Y,'single') ) Y=double(Y); end;
if(isa(K,'single')) eps=1e-7; else eps=1e-16; end;

wb=opts.alphab; if ( isempty(wb) ) wb=zeros(N+1,1); end 

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
err  = 1-Y'.*(wK+wb(end)); svs=err>=0 & Y'~=0;
% pre-condinationed gradient, 
Yerr = wghtY'.*err; % wght.*Y.*(1-Y.*f)
% K^-1*dJdw = K^-1(2 C(1)Kw - 2 K I_sv(Y-f)) = 2*(C(1)w - Isv (Y-f) )
% N.B. could include additional diag(K) pre-cond here if it would help
MdJ  = [(2*C(1)*wb(1:end-1) - 2*(Yerr.*svs)')./diag(K);...
        -2*sum(Yerr(svs))/N];
dJ   = [K*(MdJ(1:end-1).*diag(K))+opts.ridge.*MdJ(1:end-1).*diag(K); ...
        N*MdJ(end)];
if ( opts.nobias ) MdJ(end)=0; dJ(end)=0; end;
Mr   =-MdJ;
d    = Mr;
ddJ  =-d'*dJ;
r2   = ddJ;
r02  = r2;

Ed   = (Y(svs)'.*Yerr(svs))*err(svs)';  
Ew   = wK*wb(1:end-1);
J    = Ed + C(1)*Ew; % SVM objective

step=opts.step;
if( step<=0 ) step=abs(J/ddJ); end  % init step assuming opt is at 0
step=abs(step); tstep=step;

neval=1;
if(opts.verb>0)   % debug code      
   fprintf(['%3d) %3d x=[%8g,%8g,.] J=%5g |dJ|=%8g\n'],0,neval,wb(1),wb(2),J,r2);
end

% alt pre-cond non-lin CG iteration
fin1=0;
for i=1:opts.maxIter;

   oJ= J; or2=r2; oMr=Mr; % record info about prev result we need

   %---------------------------------------------------------------------
   % Secant method for the root search.
   if ( opts.verb > 1 )
      fprintf('.%d %g=%g @ %g\n',0,0,ddJ,J);
      if ( opts.verb>2 ) 
         hold off;plot(0,ddJ,'r*');hold on;text(0,double(ddJ),'0');
      end
   end;
   step=max(tstep,1e-6/ddJ); % prev step size is first guess!
   oddJ=ddJ; % one step before is same as current
   wb = wb + step*d;
   Kd = [K*d(1:end-1)+opts.ridge'.*d(1:end-1);d(end)];% cache so don't comp dJ
   for j=1:50;
      neval = neval+1;      
      ooddJ=oddJ; oddJ=ddJ; % prev and 1 before grad values
      
      % Eval the gradient at this point.  N.B. only gradient needed for secant
      wK   = wb(1:end-1)'*K + opts.ridge'.*wb(1:end-1)'; % include ridge
      err  = 1-Y'.*(wK+wb(end)); svs=err>0 & Y'~=0;
      Yerr = wghtY'.*err; 
      MdJ  = [(2*C(1)*wb(1:end-1) - 2*(Yerr.*svs)'); ...
              -2*sum(Yerr(svs))];
      if ( opts.nobias ) MdJ(end)=0; end;
      ddJ  =-Kd'*MdJ;  % gradient along the line      
      
      if ( opts.verb > 1 )
         Ed     = err(svs)*err(svs)';  Ew     = wK*wb(1:end-1);
         J      = Ed + C(1)*Ew; % SVM objective
         fprintf('.%d %g=%g @ %g \n',j,tstep,ddJ,J); 
         if ( opts.verb > 2 )
            plot(tstep,ddJ,'*'); text(double(tstep),double(ddJ),num2str(j));
         end;
      end;
      
      % convergence test, and numerical res test
      if ( abs(ddJ) <= opts.tol || abs(ddJ*step)<eps ) break; end; 

      % now compute the new step size
      % backeting check, so it always decreases
      if ( ooddJ*oddJ < 0 && oddJ*ddJ > 0 ...     % ooddJ still brackets
           && abs(step*oddJ/(oddJ-ddJ)) > abs(ostep) ) % would jump outside 
         step = ostep + step; % make as if we jumped here directly.
         oddJ = ooddJ;
      end
      ostep = step;
      % *RELATIVE* secant step size
      nstep = ddJ/(oddJ-ddJ);
      nstep=sign(nstep)*max(opts.minStep,min(abs(nstep),opts.maxStep));
      step  = step * nstep;
      tstep = tstep + step;            % total step size
      
      % move to the new point
      wb    = wb + step*d ;
   end
   if ( opts.verb > 1 ) fprintf('\n'); end;
   % compute the other bits needed for CG iteration
   dJ   = [K*(MdJ(1:end-1))+opts.ridge.*(MdJ(1:end-1)); ...
           MdJ(end)];
   if ( opts.nobias ) dJ(end)=0; end;
   Mr    =-MdJ;
   r2    =abs(Mr'*dJ);
      
   delta = max((Mr-oMr)'*(-dJ)/or2,0); % Polak-Ribier
   %delta = max(r2/or2,0); % Fletcher-Reeves
   d     = Mr+delta*d;    % conj grad direction
   ddJ   =-d'*dJ;         % new search dir grad.
   if( ddJ <= 0 )         % non-descent dir switch to steepest
      if ( opts.verb >= 1 ) fprintf('non-descent dir\n'); end;      
      d=Mr; ddJ=-d'*dJ; 
   end; 

   if ( opts.weightDecay > 0 ) % decay term to 0 non-svs faster
      pts=~svs' & wb(1:end-1).*d(1:end-1) < 0;
      %wb(pts)=wb(pts)*opts.weightDecay; 
      d(pts) =d(pts)*opts.weightDecay;
   end
   
   % compute the function evaluation
   Ed     = (Y(svs)'.*Yerr(svs))*err(svs)';  
   Ew     = wK*wb(1:end-1);
   J      = Ed + C(1)*Ew; % SVM objective
   if(opts.verb>0)   % debug code      
      fprintf(['%3d) %3d x=[%8g,%8g,.] J=%5g |dJ|=%8g\n'],...
              i,neval,wb(1),wb(2),J,r2);
   end   

   if ( r2<=opts.tol || r2<=r02*opts.tol0 || ...
        neval > opts.maxEval || abs(oJ-J) <= opts.objTol ...
        || oJ-J < -1e-2) % term tests
      break;
   end; 

end;
if opts.verb >= 0
    fprintf(['%3d) %3d x=[%8g,%8g,.] J=%5g |dJ|=%8g\n'],...
        i,neval,wb(1),wb(2),J,r2);
end

% compute final decision values.
f = wb(1:end-1)'*K + wb(end);

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
return

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

