function [apcs,apcss]=gp_avpredcomp(gp, x, y, varargin)
%GP_AVPREDCOMP  Average predictive comparison for Gaussian process model
%
%  Description
%    APCS=GP_AVPREDCOMP(GP, X, Y, OPTIONS) Takes a Gaussian process
%    structure GP together with a matrix X of training inputs and
%    vector Y of training targets, and returns average predictive
%    comparison (APC) estimates for each input in a structure APCS. 
%    APCS contains following fields
%      ps    - the probability of knowing the sign of the APC
%              in the latent outcome for each input variable.
%      fs    - the samples from the APC in the latent outcome for each
%              input variable
%      fsa   - the samples from the absolute APC in the latent outcome 
%              for each input variable
%      fsrms - the samples from the root mean squared APC in the latent
%              outcome for each input variable
%      ys    - the samples from the APC in the target outcome for each
%              input variable
%      ysa   - the samples from the absolute APC in the target outcome 
%              for each input variable
%      ysrms - the samples from the root mean squared APC in the target
%              outcome for each input variable
%
%    [APCS,APCSS]=GP_AVPREDCOMP(GP, X, Y, OPTIONS) returns also APCSS
%    which contains APCS components for each data point. These can
%    be used to form conditional average predictive comparisons (CAPC).
%    APCSS contains following fields
%      numfs    - the samples from the numerator of APC in the latent
%                 outcome for each input variable
%      numfsa   - the samples from the numerator of absolute APC in
%                 the latent outcome for each input variable
%      numfsrms - the samples from the numerator of RMS APC in
%                 the latent outcome for each input variable
%      numys    - the samples from the numerator of APC in the latent
%                 outcome for each input variable
%      numysa   - the samples from the numerator of absolute APC in
%                 the latent outcome for each input variable
%      numysrms - the samples from the numerator of RMS APC in
%                 the latent outcome for each input variable
%      dens     - the samples from the denominator of APC in the latent
%                 outcome for each input variable
%      densa    - the samples from the denominator of absolute APC in
%                 the latent outcome for each input variable
%      densrms  - the samples from the denominator of RMS APC in
%                 the latent outcome for each input variable
%
%    OPTIONS is optional parameter-value pair
%      z         - optional observed quantity in triplet (x_i,y_i,z_i)
%                  Some likelihoods may use this. For example, in
%                  case of Poisson likelihood we have z_i=E_i, that
%                  is, expected value for ith case.
%      nsamp     - determines the number of samples used (default=500).
%      deltadist - indicator vector telling which component sets
%                  are handled using the delta distance (0 if x=x',
%                  and 1 otherwise). Default is found by examining
%                  the covariance and metric functions used.
%
%  See also
%    GP_PRED

% Copyright (c) 2011      Jaakko Riihimäki
% Copyright (c) 2011      Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.FunctionName = 'GP_AVPREDCOMP';
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('predcf', [], @(x) isempty(x) || ...
                 isvector(x) && isreal(x) && all(isfinite(x)&x>0))
ip.addParamValue('tstind', [], @(x) isempty(x) || iscell(x) ||...
                 (isvector(x) && isreal(x) && all(isfinite(x)&x>0)))
ip.addParamValue('nsamp', 500, @(x) isreal(x) && isscalar(x))
ip.addParamValue('deltadist',[], @(x) isvector(x));

ip.parse(gp, x, y, varargin{:});
options=struct();
options.predcf=ip.Results.predcf;
options.tstind=ip.Results.tstind;
z=isempty(ip.Results.z);
if ~isempty(z)
  options.z=ip.Results.z;
end
nsamp=ip.Results.nsamp;
deltadist = logical(ip.Results.deltadist);

[n, nin]=size(x);
if isempty(deltadist)
  deltadist=false(1,nin);
  deltadist(gp_finddeltadist(gp))=true;
end

ps=zeros(1,nin);
fs=zeros(nsamp,nin);
fsa=zeros(nsamp,nin);
fsrms=zeros(nsamp,nin);
if nargout>1
  numfs=zeros(n,nsamp,nin);
  numfsa=zeros(n,nsamp,nin);
  numfsrms=zeros(n,nsamp,nin);
  dens=zeros(n,nsamp,nin);
  densa=zeros(n,nsamp,nin);
  densrms=zeros(n,nsamp,nin);
end
ys=zeros(nsamp,nin);
ysa=zeros(nsamp,nin);
ysrms=zeros(nsamp,nin);
if nargout>1
  numys=zeros(n,nsamp,nin);
  numysa=zeros(n,nsamp,nin);
  numysrms=zeros(n,nsamp,nin);
end

% covariance is used for Mahalanobis weighted distanec
covx=cov(x);
% handle categorical variables
covx(deltadist,:)=0;
covx(:,deltadist)=0;
for i1=find(deltadist)
  covx(i1,i1)=1;
end

prevstream=setrandstream();

% loop through the input variables
for k1=1:nin
  fprintf('k1=%d\n',k1)
  %- Compute the weight matrix based on Mahalanobis distances:
  x_=x; x_(:,k1)=[];
  covx_=covx; covx_(:,k1)=[]; covx_(k1,:)=[];
  deltadist_=deltadist; deltadist_(k1)=[];
  % weight matrix:
  W=zeros(n);
  for i1=1:n
    x_diff=zeros(nin-1,n-i1);
    x_diff(~deltadist_,:)=bsxfun(@minus,x_(i1,~deltadist_),x_((i1+1):n,~deltadist_))';
    x_diff(deltadist_,:)=double(bsxfun(@ne,x_(i1,deltadist_),x_((i1+1):n,deltadist_))');
    W(i1,(i1+1):n)=1./(1+sum(x_diff.*(covx_\x_diff)));
  end
  W=W+W'+eye(n);
  
  seed=round(rand*10e8);
  
  numf=zeros(1,nsamp);
  numfa=zeros(1,nsamp);
  numfrms=zeros(1,nsamp);
  numy=zeros(1,nsamp); 
  numya=zeros(1,nsamp);
  numyrms=zeros(1,nsamp);
  den=0;
  dena=0;
  for i1=1:n
    % inputs of interest
    ui=x(i1, k1);
    ujs=x(:, k1);
    
    % replicate same values for other inputs
    xrep=repmat(x(i1,:),n,1); xrep(:,k1)=ujs;
    
    if deltadist(k1)
      Udiff=double(ujs~=ui);
    else
      Udiff=ujs-ui;
    end
    Udiffa=abs(Udiff);
    Usign=sign(Udiff);

    % draw random samples from the posterior
    setrandstream(seed);
    fr = gp_rnd(gp, x, y, xrep, 'nsamp', nsamp, options);
    
    % average change in input
    deni=sum(W(:,i1).*Udiff.*Usign);
    denai=sum(W(:,i1).*Udiffa);
    den=den+deni;
    dena=dena+denai;
    
    % average change in latent outcome
    b=bsxfun(@minus,fr,fr(i1,:));
    numfi=sum(bsxfun(@times,W(:,i1).*Usign,b));
    numfai=sum(bsxfun(@times,W(:,i1),abs(b)));
    numfrmsi=sum(bsxfun(@times,W(:,i1),b.^2));
    numf=numf+numfi;
    numfa=numfa+numfai;
    numfrms=numfrms+numfrmsi;

    if nargout>1
      numfs(i1,:,k1)=numfi;
      numsa(i1,:,k1)=numfai;
      numfsrms(i1,:,k1)=numfrmsi;
      dens(i1,:,k1)=deni;
      densa(i1,:,k1)=denai;
      densrms(i1,:,k1)=denai;
    end
    
    % compute latent values through the inverse link function
    if isfield(gp.lik.fh, 'invlink')
      ilfr = gp.lik.fh.invlink(gp.lik, fr, repmat(z,1,nsamp));
      % average change in outcome
      b=bsxfun(@minus,ilfr,ilfr(i1,:));
      numyi=sum(bsxfun(@times,W(:,i1).*Usign,b));
      numyai=sum(bsxfun(@times,W(:,i1),abs(b)));
      numyrmsi=sum(bsxfun(@times,W(:,i1),b.^2));
      numy=numy+numyi;
      numya=numya+numyai;
      numyrms=numyrms+numyrmsi;
      
      if nargout>1
        numys(i1,:,k1)=numyi;
        numysa(i1,:,k1)=numyai;
        numysrms(i1,:,k1)=numyrmsi;
      end
    end
    
  end
  
  % outcome is the latent function
  fs(:,k1)=numf./den;
  fsa(:,k1)=numfa./dena;
  fsrms(:,k1)=sqrt(numfrms./dena);
  
  if isfield(gp.lik.fh, 'invlink')
    % outcome is computed through the inverse link function
    ys(:,k1)=numy./den;
    ysa(:,k1)=numya./dena;
    ysrms(:,k1)=sqrt(numyrms./dena);
  end

  % probability of knowing the sign of the change in
  % latent function 
  ps(1,k1)=mean(numf./den>0);
  if ps(1,k1)<0.5
    ps(1,k1)=1-ps(1,k1);
  end
  
end

apcs.ps=ps;
apcs.fs=fs;
apcs.fsa=fsa;
apcs.fsrms=fsrms;
if isfield(gp.lik.fh, 'invlink')
  apcs.ys=ys;
  apcs.ysa=ysa;
  apcs.ysrms=ysrms;
end

if nargout>1
  apcss.numfs=numfs;
  apcss.numfsa=numfsa;
  apcss.numfsrms=numfsrms;
  apcss.dens=dens;
  apcss.densa=densa;
  apcss.densrms=densrms;
  if isfield(gp.lik.fh, 'invlink')
    apcss.numys=numys;
    apcss.numysa=numysa;
    apcss.numysrms=numysrms;
  end
end

setrandstream(prevstream);

end

function deltadist = gp_finddeltadist(cf)
% FINDDELTADIST - Find which covariates are using delta distance
%   
deltadist=[];
if ~iscell(cf) && isfield(cf,'cf')
  deltadist=union(deltadist,gp_finddeltadist(cf.cf));
else
  for cfi=1:numel(cf)
    if isfield(cf{cfi},'cf')
      deltadist=union(deltadist,gp_finddeltadist(cf{cfi}.cf));
    else
      if isfield(cf{cfi},'metric')
        if isfield(cf{cfi}.metric,'deltadist')
          deltadist=union(deltadist,cf{cfi}.metric.deltadist);
        end
      elseif ismember(cf{cfi}.type,{'gpcf_cat' 'gpcf_mask'}) && ...
          isfield(cf{cfi},'selectedVariables')
        deltadist=union(deltadist,cf{cfi}.selectedVariables);
      end
    end
  end
end
end
