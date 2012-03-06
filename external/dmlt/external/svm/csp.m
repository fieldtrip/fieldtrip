function [sf,d,Sigmai,Sigmac,SigmaAll]=csp(X,Y,dim,cent,ridge,singThresh)
% Generate spatial filters using CSP
%
% [sf,d,Sigmai,Sigmac,SigmaAll]=csp(X,Y,[dim]);
% N.B. if inputs are singular then d will contain 0 eigenvalues & sf==0
% Inputs:
%  X     -- n-d data matrix, e.g. [nCh x nSamp x nTrials] data set, OR
%           [nCh x nCh x nTrials] set of *trial* covariance matrices, OR
%           [nCh x nCh x nClass ] set of *class* covariance matrices
%  Y     -- [nTrials x 1] set of trial labels, with nClass unique labels, OR
%           [nTrials x nClass] set of +/-1 (&0) trial lables per class, OR
%           [nClass  x 1] set of class labels when X=[nCh x nCh x nClass] (1:size(X,dim))
%           N.B. in all cases a label of 0 indicates ignored trial
%  dim   -- [1 x 2] dimension of X which contains the trials, and
%           (optionally) the the one which contains the channels.  If
%           channel dim not given the next available dim is used. ([-1 1])
%  cent  -- [bool] center the data (0)
%  ridge -- [float] size of ridge (as fraction of mean eigenvalue) to add for numerical stability (1e-7)
%  singThresh -- [float] threshold to detect singular values in inputs (1e-3)
% Outputs:
%  sf    -- [nCh x nCh x nClass] sets of 1-vs-rest spatial *filters*
%           sorted in order of increasing eigenvalue.
%           N.B. sf is normalised such that: mean_i sf'*cov(X_i)*sf = I
%           N.B. to obtain spatial *patterns* just use, sp = Sigma*sf ;
%  d     -- [nCh x nClass] spatial filter eigen values
%  Sigmai-- [nCh x nCh x nTrials] set of *trial* covariance matrices
%  Sigmac-- [nCh x nCh x nClass]  set of *class* covariance matrices
%  SigmaAll -- [nCh x nCh] all (non excluded) data covariance matrix
%
if ( nargin < 3 || isempty(dim) ) dim=[-1 1]; end;
if ( numel(dim) < 2 ) if ( dim(1)==1 ) dim(2)=2; else dim(2)=1; end; end
dim(dim<0)=ndims(X)+dim(dim<0)+1; % convert negative dims
if ( nargin < 4 || isempty(cent) ) cent=0; end;
if ( nargin < 5 || isempty(ridge) ) 
   if ( isequal(class(X),'single') ) ridge=1e-7; else ridge=0; end;
end;
if ( nargin < 6 || isempty(singThresh) ) singThresh=1e-3; end

nCh = size(X,dim(2)); N=size(X,dim(1)); nSamp=prod(size(X))./nCh./N;

% compute the per-trial covariances
if ( ~isequal(dim,[3 1]) || ndims(X)>3 || nCh ~= size(X,2) )         
   idx1=-(1:ndims(X)); idx2=-(1:ndims(X)); % sum out everything but ch, trials
   idx1(dim(1))=3;     idx2(dim(1))=3;     % linear over trial dimension
   idx1(dim(2))=1;     idx2(dim(2))=2;     % Outer product over ch dimension
   Sigmai = tprod(X,idx1,[],idx2,'n');
   if ( cent ) % center the co-variances, N.B. tprod to comp means for mem
      error('Unsupported -- numerically unsound, center before instead');
%       sizeX=size(X); muSz=sizeX; muSz(dim)=1; mu=ones(muSz,class(X));
%       idx2(dim)=0; mu=tprod(X,idx1,mu,idx2); % nCh x 1 x nTr
%       % subtract the means
%       Sigmai = Sigmai - tprod(mu,[1 0 3],[],[2 0 3])/prod(muSz); 
   end
%    Fallback code
%    if(dim(1)==3)     for i=1:size(X,3); Sigmai(:,:,i)=X(:,:,i)*X(:,:,i)'; end
%    elseif(dim(1)==1) for i=1:size(X,1); Sigmai(:,:,i)=shiftdim(X(i,:,:))*shiftdim(X(i,:,:))'; end
%    end
else
   Sigmai = X;
end
% N.B. Sigmai *must* be [nCh x nCh x N]

if ( ndims(Y)==2 && min(size(Y))==1 && ~(all(Y(:)==-1 | Y(:)==0 | Y(:)==1)) ) 
  oY=Y;
  Y=lab2ind(Y,[],[],[],0); 
end;
nClass=size(Y,2);
if ( nClass==2 ) nClass=1; end; % only 1 for binary problems

allY0 = all(Y==0,2); % trials which have label 0 in all sub-prob
SigmaAll = sum(double(Sigmai(:,:,~allY0)),3); % sum all non-0 labeled trials
sf    = zeros([nCh,nCh,nClass],class(X)); d=zeros(nCh,nClass,class(X));
for c=1:nClass; % generate sf's for each sub-problem
   Sigmac(:,:,c) = sum(double(Sigmai(:,:,Y(:,c)>0)),3); % +class covariance
   if ( isequal(Y(:,c)==0,allY0) ) % rest covariance, i.e. excludes 0 class
      Sigma=SigmaAll; % can use sum of everything
   else
      Sigma=sum(Sigmai(:,:,Y(:,c)~=0),3); % rest for this class
   end
   Sigma=double(Sigma);
   
   % solve the generalised eigenvalue problem, 
   if ( ridge>0 ) % Add ridge if wanted to help with numeric issues in the inv
      Sigmac(:,:,c)=Sigmac(:,:,c)+eye(size(Sigma))*ridge*mean(diag(Sigmac(:,:,c))); 
      Sigma        =Sigma        +eye(size(Sigma))*ridge*mean(diag(Sigma));
   end
   % N.B. use double to avoid rounding issues with the inv(Sigma) bit
   [W D]=eig(Sigmac(:,:,c),Sigma);D=diag(D);
   [dc,di]=sort(D,'descend');   W=W(:,di); % order in decreasing eigenvalue

   % Check for and correct for singular inputs
   % singular if eigval out of range
   nf = sum(W.*(Sigma*W),1)'; % eig-value for this direction in full cov
   si= dc>1-singThresh | dc<0+singThresh | imag(dc)~=0 | isnan(dc) | nf<1e-4*sum(abs(nf)); 
   if ( sum(si)>0 ) % remove the singular eigen values & effect on other sf's

      % Identify singular directions which leak redundant information into
      % the other eigenvectors, i.e. are mapped to 0 by both Sigmac and Sigma
      % N.B. if the numerics are OK this is probably uncessary!
      Na  = sum((double(Sigmac(:,:,c))*W(:,si)).^2)./sum(W(:,si).^2); 
      ssi=find(si);ssi=ssi(abs(Na)<singThresh & imag(Na)==0);%ssi=dc>1-singThresh|dc<0+singThresh;
      if ( ~isempty(ssi) ) % remove anything in this dir in other eigenvectors
         % Compute the projection of the rest onto the singular direction(s)
         Pssi    = repop(W(:,ssi)'*W(:,~si),'./',sum(W(:,ssi).*W(:,ssi),1)');
         W(:,~si)= W(:,~si) - W(:,ssi)*Pssi; %remove this singular contribution
      end

      W=W(:,~si); dc=dc(~si); nf=nf(~si);% discard singular components   
   end
   
   %Normalise, so that diag(W'*Sigma*W)=N, i.e.mean_i W'*(X_i*X_i')*W/nSamp = 1
   % i.e. so that the resulting features have unit variance (and are
   % approx white?)
   W = repop(W,'*',nf'.^-.5)*sqrt(sum(Y(:,c)~=0)*nSamp); 

   % Save the normalised filters & eigenvalues
   sf(:,1:size(W,2),c)= W;  
   d(1:size(W,2),c)   = dc;
end
% Compute last class covariance if wanted
if ( nClass==1 & nargout>3 ) Sigmac(:,:,2)=sum(double(Sigmai(:,:,Y(:,1)<0)),3);end;
return;

%-----------------------------------------------------------------------------
function []=testCase()
nCh = 64; nSamp = 100; N=300;
X=randn(nCh,nSamp,N);
Y=sign(randn(N,1));
[sf,d,Sigmai,Sigmac]=jf_csp(X,Y);

[sf,d,Sigmai,Sigmac]=jf_csp(X,Y,1); % with data centering

[sf2,d2]=jf_csp(Sigmac,[-1 1]); 
[sf3,d3]=csp(Sigmac);
mimage(sf,sf2,'diff',1,'clim','limits')