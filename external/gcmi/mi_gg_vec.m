function I = mi_gg_vec(x, y, biascorrect, demeaned)
% MI_GG Vectorized MI calculation between multiple Gaussian x variables and a 
%   common Gaussian y variable in bits
%   I = mi_gg_vec(x,y) returns the MI between two (possibly multidimensional)
%   Gassian variables, x and y, with bias correction.
%   size(x) = [Ntrl Nvec Ndim]
%   so each output I(i) = mi_gg_vec(squeeze(x(:,i,:)), y);
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
if ndims(x)>3 
  error('mi_gg_vec: x input array should be 3d')
end
if ndims(y)~=2
  error('mi_gg_vec: y input should be 2d')
end

Ntrl = size(x,1);
Nvec = size(x,2);
Nvarx = size(x,3);
Nvary = size(y,2);
Nvarxy = Nvarx + Nvary;

if size(y,1) ~= Ntrl
    error('mi_gg_vec: number of trials do not match')
end

% default option values
if nargin<3
    biascorrect = true;
end
if nargin<4
    demeaned = false;
end

% demean data if required
if ~demeaned
    x = bsxfun(@minus,x,sum(x,1)/Ntrl);
    y = bsxfun(@minus,y,sum(y,1)/Ntrl);
end

Cx = zeros(Nvec, Nvarx, Nvarx);
Cxy = zeros(Nvec, Nvarxy, Nvarxy);

Cy = y'*y / (Ntrl - 1);

% Cx and Cx part of Cxy
for vi1=1:Nvarx
  x1 = x(:,:,vi1);
  thsV = sum(x1.^2);
  Cx(:,vi1,vi1) = thsV;
  Cxy(:,vi1,vi1) = thsV;
  
  for vi2=(vi1+1):Nvarx
    x2 = x(:,:,vi2);
    thsC = sum(x1.*x2);
    Cx(:,vi1,vi2) = thsC;
    Cx(:,vi2,vi1) = thsC;
    Cxy(:,vi1,vi2) = thsC;
    Cxy(:,vi2,vi1) = thsC;
  end
end

Cx = Cx / (Ntrl-1);

% Cxy part of Cxy
for vi1=1:Nvarx
  x1 = x(:,:,vi1);
  for vi2=1:Nvary
    y1 = y(:,vi2);
    thsC = sum(bsxfun(@times,x1,y1));
    Cxy(:,Nvarx+vi2,vi1) = thsC;
    Cxy(:,vi1,Nvarx+vi2) = thsC;
  end
end

Cxy = Cxy ./ (Ntrl-1);

% Cy part of Cxy
for vi1=1:Nvary
  Cxy(:,Nvarx+vi1,Nvarx+vi1) = Cy(vi1,vi1);
  for vi2=(vi1+1):Nvary
    Cxy(:,Nvarx+vi1,Nvarx+vi2) = Cy(vi1,vi2);
    Cxy(:,Nvarx+vi2,Nvarx+vi1) = Cy(vi1,vi2);
  end
end


% entropies in nats
% normalisations cancel for information
chCy = chol(Cy);
HY = sum(log(diag(chCy))); % + 0.5*Nvary*log(2*pi*exp(1));

chCx = vecchol(Cx);
chCxy = vecchol(Cxy);
HX = zeros(Nvec,1);
HXY = zeros(Nvec,1);
for vi=1:Nvarx
%   HX = HX + shiftdim(log(Cx(:,vi,vi)));
  HX = HX + log(chCx(:,vi,vi));
end
for vi=1:Nvarxy
%   HXY = HXY + shiftdim(log(Cxy(:,vi,vi)));
  HXY = HXY + log(chCxy(:,vi,vi));
end

ln2 = log(2);
if biascorrect
    psiterms = psi((Ntrl - (1:Nvarxy))/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
%     HX = (HX - Nvarx*dterm - sum(psiterms(1:Nvarx)));
%     HY = (HY - Nvary*dterm - sum(psiterms(1:Nvary)));
%     HXY = (HXY - Nvarxy*dterm - sum(psiterms));
    HXbias = Nvarx*dterm + sum(psiterms(1:Nvarx));
    HYbias = Nvary*dterm + sum(psiterms(1:Nvary));
    HXYbias = Nvarxy*dterm + sum(psiterms);
    Ibias = HXbias + HYbias - HXYbias;
else
    Ibias = 0;
end

% convert to bits
% I = (HX + HY - HXY) / ln2;
I = (bsxfun(@plus,HX-HXY, HY) - Ibias) ./ ln2;

