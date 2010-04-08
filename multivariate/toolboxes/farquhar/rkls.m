function [wb,f,J]=rkls(K,Y,C,varargin)
% Regularised Kernel Least Squares Classifier
%
% [alphab,J]=rkls(K,Y,C,varargin)
% Simple regularised kernel least squares classifier.
%
% J = C(1)*w'*K*w + sum_i (y_i - (w'*K_i + b)).^2
%
% Inputs:
% K       - [N x N] kernel matrix
% Y       - [N x L] matrix of +1, -1 labels, for an L-class problem,
%           N.B. points with label 0 are ignored
% C       - [1 x 1] regularisation parameter
%
% Outputs:
% alphab  - [(N+1) x 1] matrix of kernel weights and bias [alpha;b]
% f       - [N x 1] The decision value for all the inputs
% J       - the final objective value
%
% Options:
%
%  wght    - point weights [Nx1] vector of label accuracy probabilities
%            [2x1] for per class weightings
%  nobias  - flag we don't want the bias computed (false)
%  verb    - verbosity level (0)
% 
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)

% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied

% Argument processing
if ( nargin < 3 ) C(1)=0; end;
opts=struct('alphab',[],'nobias',0,'verb',0);
[opts,varargin]=parseOpts(opts,varargin{:});

% Identify the training points
trnInd = (Y~=0); 
N=sum(trnInd);

% Add the bias term to the kernel -- N.B. this is a regularised bias.
if ( ~opts.nobias ) biasC = 1; else biasC = 0 ; end

% Now train on the training points.
w = [C*eye(N) + K(trnInd,trnInd) + biasC]\Y(trnInd) ;
% Extract the bias term from the wb
b = sum(w);  % N.B. should be 0 if nobias is set
% The final solution
wb = zeros(size(K,1)+1,1); % N.B. w.r.t. the full input set
wb(trnInd) = w; wb(end)=b;

% Now compute the predictions and the objective function value
Kw = K*wb(1:end-1);
f  = Kw + b;
J  = C(1)*wb(1:end-1)'*Kw + sum((Y(trnInd)-f(trnInd)).^2);
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

%-------------------------------------------------------------------------
function testCase()

[X,Y]=mkMultiClassTst([-1 0; 1 0; .2 .5],[400 400 50],[.3 .3; .3 .3; .2 .2],[],[-1 1 1]);[dim,N]=size(X);

K=X'*X; % N.B. add 1 to give implicit bias term
fInds=gennFold(Y,10,'perm',1); trnInd=any(fInds(:,1:9),2); tstInd=fInds(:,10);

[alphab,f,J]=rkls(K,Y.*double(trnInd),1);
dv=K*alphab(1:end-1)+alphab(end);
dv2conf(Y(tstInd),dv(tstInd))

% for linear kernel
plotLinDecisFn(X,Y,X*alphab(1:end-1),alphab(end),alpha);

