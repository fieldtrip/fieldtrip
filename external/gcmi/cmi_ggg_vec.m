function I = cmi_ggg_vec(x, y, z, biascorrect, demeaned)
% CMI_GGG Conditional Mutual information (CMI) between two Gaussian variables
%        conditioned on a third
%
%   I = cmi_ggg(x,y,z) returns the CMI between two (possibly multidimensional)
%   Gassian variables, x and y, conditioned on a third, z, with bias correction.
%   If x / y / z are multivariate rows must correspond to samples, columns
%   to dimensions/variables. (Samples first axis) 
%
%   biascorrect : true / false option (default true) which specifies whether
%   bias correction should be applied to the esimtated MI.
%   demeaned : false / true option (default false) which specifies whether the
%   input data already has zero mean (true if it has been copula-normalized)

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if isvector(y)
    y = y(:);
end
if isvector(z)
    z = z(:);
end

[Ntrlx,Nsgnx,Nvarx] = size(x);
[Ntrly,Nsgny,Nvary] = size(y);
[Ntrlz,Nsgnz,Nvarz] = size(z);

% check here if Nsgnz==Nsgnx
if Nsgnx==Nsgnz
  % ok
elseif Nsgnx==1 && Nsgnz>1
  x = repmat(x, [1 Nsgnz]);
  Nsgnx = size(x,2);
elseif Nsgnz==1 && Nsgnx>1
  z = repmat(z, [1 Nsgnx]);
  Nsgnz = size(z,2);
else
  error('unsupported dimensionality of input');
end

if (Ntrlx~=Ntrly) || (Ntrlx~=Ntrlz)
    error('cmi_ggg: number of trials do not match')
end

% default option values
if nargin<4
    biascorrect = true;
end
if nargin<5
    demeaned = false;
end

% demean data if required
if ~demeaned
    x = bsxfun(@minus,x,sum(x,1)/Ntrlx);
    y = bsxfun(@minus,y,sum(y,1)/Ntrly);
    z = bsxfun(@minus,z,sum(z,1)/Ntrlz);
end

Nvaryz  = Nvary + Nvarz;
Nvarxyz = Nvarx + Nvaryz;
Nvarxz  = Nvarx + Nvarz;

xidx = 1:Nvarx;
zidx = Nvarx+Nvary+(1:Nvarz);
idx  = (Nvarx + 1):Nvarxyz;

% joint variable
Cxyz = zeros(Nsgnx,Nvarxyz,Nvarxyz);

for k = 1:Nsgnx
  xyz = [reshape(x(:,k,:),[Ntrlx Nvarx]) reshape(y, [Ntrly Nvary*Nsgny]) reshape(z(:,k,:), [Ntrlz Nvarz])];
  Cxyz(k,:,:) = (xyz'*xyz) / (Ntrlx - 1);
end

% submatrices of joint covariance
Cz  = Cxyz(:,zidx, zidx);
Cyz = Cxyz(:,idx,  idx);


Cxz = zeros(Nsgnx, Nvarxz, Nvarxz);
Cxz(:,xidx,xidx) = Cxyz(:,xidx,xidx);
zidxxz = (Nvarx+1):Nvarxz;
Cxz(:,xidx,zidxxz) = Cxyz(:,xidx,zidx);
Cxz(:,zidxxz,xidx) = Cxyz(:,zidx,xidx);
Cxz(:,zidxxz,zidxxz) = Cxyz(:,zidx,zidx);

chCz   = vecchol(Cz);
chCxz  = vecchol(Cxz);
chCyz  = real(vecchol(Cyz));
chCxyz = real(vecchol(Cxyz));

% entropies in nats
% normalisations cancel for cmi
HZ   = sum(log(vecdiag(chCz)),  2); % + 0.5*Nvarz*log(2*pi*exp(1));
HXZ  = sum(log(vecdiag(chCxz)), 2); % + 0.5*(Nvarx+Nvarz)*log(2*pi*exp(1));
HYZ  = sum(log(vecdiag(chCyz)), 2); % + 0.5*(Nvary+Nvarz)*log(2*pi*exp(1));
HXYZ = sum(log(vecdiag(chCxyz)),2); % + 0.5*(Nvarx+Nvary+Nvarz)*log(2*pi*exp(1));

ln2 = log(2);
if biascorrect
    psiterms = psi((Ntrlx - (1:Nvarxyz))/2) / 2;
    dterm = (ln2 - log(Ntrlx-1)) / 2;
    HZ   = (HZ   - Nvarz*dterm   - sum(psiterms(1:Nvarz)));
    HXZ  = (HXZ  - Nvarxz*dterm  - sum(psiterms(1:Nvarxz)));
    HYZ  = (HYZ  - Nvaryz*dterm  - sum(psiterms(1:Nvaryz)));
    HXYZ = (HXYZ - Nvarxyz*dterm - sum(psiterms));
end

% convert to bits
I = (HXZ + HYZ - HXYZ - HZ) / ln2;

function out = vecdiag(in)

n = size(in,2);
out = zeros(size(in,1),n);
for k = 1:n
  out(:,k) = in(:,k,k);
end
