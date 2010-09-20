function [beta,yhat,obfun]=elasticlr(X,y,options)
% Elastic Net Logistic Regression for two class problem using
% coordinate descent algorithm.
% The implementation is based on the paper by Friedman, Hastie & Tibshirani,
% "Regularization paths for generalized linear models via coordinate descent",
% Journal of Statistical software, Vol.33(1), pp 1-22
%
% X: ntrial x nfeatures (assumed to be standardized!!!)
% y: ntrial x 1
% beta: (nfeatures+1) x 1
% options.L1=0.1;      %L1 regularization parameter
% options.L2=0.1;      %L2 regularization parameter
% options.eta=0.1;     %stepsize for the Newton's update
% options.HessApprox=1;%use Hessian approximation
% options.tol=1e-6;    %tolerance of changes in Beta
% options.maxiter=5e3; %maximum number of iterations
% options.verbose=1;   %show output
%
% bias term not included by default!
%
% all terms (including possibly bias) are regularized
%
% to do: allow ridge penalty to induce spatial smoothness
%
% copyright 2010, Ali Bahramisharif

  if nargin<2
    error('data and labels are required')
  end
  if ~exist('options','var')
    options=[];
  end
  if ~isfield(options,'L1')
    options.L1=0.1;
  end
  if ~isfield(options,'L2')
    options.L2=0.1;
  end
  if ~isfield(options,'eta')
    options.eta=0.01;
  end
  if ~isfield(options,'maxiter')
    options.maxiter=1e4;
  end
  if ~isfield(options,'verbose')
    options.verbose=1;
  end
  if ~isfield(options,'HessApprox')
    options.HessApprox=1;
  end
  if ~isfield(options,'tol') || isempty(options.tol)
    options.tol=1e-6;
  end
  
  if max(y)>2
    error('this function only works on binary classes')
  end
  
  if length(size(X))>2
    X=reshape(X,size(X,1),[]);
  end
  ntrials=size(X,1);
  nfeatures=size(X,2);
  
  beta=zeros(size(X,2),1); % initialization with all betas to be 0
  
  newbeta=beta;
  converge=0;
  niter=0;
  XX=X.*X;
  
  if options.HessApprox
    w=0.25*ones(ntrials,1); %approximation
    wXX=w'*XX;
    wX=0.25*X;
  end
  y=y-min(y);
  options.obfun=0;
  
  if nargout>2
    obfun = nan(options.maxiter,1);
    options.obfun=1;
  end
  
  while ~converge
    
    niter=niter+1;
    Xbeta=X*beta;
    ptild=1./(1+exp(-Xbeta)); % ntrials x 1
    ptild(ptild<1e-5)=0;
    ptild(ptild>0.99999)=1;
    
    if ~options.HessApprox
      w=ptild.*(1-ptild);  % ntrials x 1
      w(w==0)=1e-5;
      r=(y-ptild)./w;
    else
      r=4*(y-ptild);
    end
    
    for j=1:nfeatures
      
      fparam=0;
      for i=1:ntrials
        if options.HessApprox
          fparam=fparam+wX(i,j)*(r(i)+X(i,j)*beta(j));
        else
          fparam=fparam+w(i)*X(i,j)*(r(i)+X(i,j)*beta(j));
        end
      end

%        if options.HessApprox
%          fparam = sum(wX(:,j).*(r+X(:,j)*beta(j)));
%        else
%          fparam = sum(w(:).*X(:,j).*(r+X(:,j)*beta(j)));
%        end

      if options.HessApprox
        newbeta(j)=(options.eta)*signz(fparam,options.L1)./(options.L2+wXX(:,j))+(1-options.eta)*beta(j);
      else
        newbeta(j)=(options.eta)*signz(fparam,options.L1)./(options.L2+w'*XX(:,j))+(1-options.eta)*beta(j);
      end
      
    end
    
    converge= (sum(abs(beta-newbeta))<options.tol || niter>options.maxiter);
    
    if options.obfun
      lexp=Xbeta;
      lexplist=find(lexp<10);
      lexp(lexplist)=log(1+exp(Xbeta(lexplist)));
      newl=-sum(mean(y.*(Xbeta))-lexp)+sum(0.5*options.L2*beta.*beta+options.alpha*abs(beta));
      obfun(niter) = newl;
    end
    
    if options.verbose
      fprintf('iteration %d of %d: %g\n',niter,options.maxiter,sum(abs(beta-newbeta)));
    end
    
    beta=newbeta;
    
  end
  
  if niter>options.maxiter && options.verbose
    fprintf('maximum number of iterations reached\n')
  end
  
  beta = beta(:);
  
  if nargout>1
    yhat=1+(ptild>0.5);
  end
  
end

function out=signz(z,gamma)

  if abs(z)>gamma
    if z>0
      out=z-gamma;
    elseif z<0
      out=z+gamma;
    else
      out=0;
    end
  else
    out=0;
  end
  
end
