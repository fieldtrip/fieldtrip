function I = mi_model_gd_vec(x, y, Ym, biascorrect, demeaned)
% MI_MODEL_GD_VEC Vectorized MI calculation between multiple Gaussian variables 
%         and a common discrete variable in bits based on ANOVA style model comparison.
%   I = mi_model_gd_vec(x,y,Ym) returns the MI between the (possibly multidimensional)
%   Gaussian variables x and the discrete variable y.
%   size(x) = [Ntrl Nvec Ndim]
%   so each output I(i) = mi_model_gd(squeeze(x(:,i,:)), y, Ym);
%
%   For 1D x this is a lower bound to the mutual information.
%   y should contain integer values in the range [0 Ym-1] (inclusive).
%
%   biascorrect : true / false option (default true) which specifies whether
%   bias correction should be applied to the esimtated MI.
%   demeaned : false / true option (default false) which specifies whether the
%   input data already has zero mean (true if it has been copula-normalized)
%   See also: MI_MODEL_GD, MI_MIXTURE_GD_VEC

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if ndims(x)>3
    error('mi_model_gd: input arrays should be 3d')
end
if isvector(y)
    y = y(:);
else
    error('mi_model_gd: only univariate discrete variable supported');
end

Ntrl = size(x,1);
Nvec = size(x,2);
Xdim = size(x,3);

if size(y,1) ~= Ntrl
    error('mi_model_gd: number of trials do not match');
end

% default option values
if nargin<4
    biascorrect = true;
end

if nargin<5
    demeaned = false;
end

% unconditional demean
if ~demeaned
    x = bsxfun(@minus,x,sum(x,1)/Ntrl);
end

I = zeros(Nvec,1);

% y = y-1;
% for vi=1:Nvec
%   I(vi) = mi_model_gd(squeeze(x(:,vi,:)),y,Ym,true,true);
% end
% 
% return

% one-hot encoding of Y
Yhot = indexed2boolean(y);

% remove class means
[Xcen class_means] = removeclassmeans(x, Yhot);

% allocate memory for class-conditional entropies and covariances
Ntrl_y = sum(Yhot);
Hcond  = zeros(Nvec,Ym);
Cm     = zeros(Nvec,Xdim,Xdim,Ym);

% allocate memory for overall entropy and covariance
Hunc = zeros(Nvec,1);
Cx   = zeros(Nvec,Xdim,Xdim);

%  data is class-demeaned, this needs to be accounted for in the
% unconditional entropies
c = diag(sqrt(Ntrl_y))*class_means.';
c = reshape(c, [Ym Nvec Xdim]);

for vi1=1:Xdim
  % all voxels for this dimension
  x1 = Xcen(:,:,vi1);
  c1 = squeeze(c(:,:,vi1));
  
  Cx(:,vi1,vi1) = sum(x1.^2)+sum(c1.^2);
  for yi=1:Ym
    tmp = x1(Yhot(:,yi),:);
    Cm(:,vi1,vi1,yi) = sum(tmp.^2);
  end
   
  for vi2=(vi1+1):Xdim
    x2  = Xcen(:,:,vi2);
    c2  = squeeze(c(:,:,vi2));
    tmp = transpose(sum(x1.*x2) + sum(c1.*c2));
    
    Cx(:,vi1,vi2) = tmp;
    Cx(:,vi2,vi1) = tmp;
    for yi=1:Ym
      tmp = transpose(sum(x1(Yhot(:,yi),:).*x2(Yhot(:,yi),:)));
      Cm(:,vi1,vi2,yi) = tmp;
      Cm(:,vi2,vi1,yi) = tmp;
    end
  end
end

Cx = Cx / (Ntrl-1);
Cx = vecchol(Cx);
for vi=1:Xdim
  Hunc = Hunc + shiftdim(log(Cx(:,vi,vi)));
end

for yi=1:Ym
  Cm(:,:,:,yi) = Cm(:,:,:,yi) / (Ntrl_y(yi) - 1);
  Cm(:,:,:,yi) = vecchol(Cm(:,:,:,yi));
  for vi=1:Xdim
    Hcond(:,yi) = Hcond(:,yi)+shiftdim(log(Cm(:,vi,vi,yi)));
  end
end

% apply bias corrections
ln2 = log(2);
if biascorrect

  vars = 1:Xdim;
  
  psiterms_unc = psi((Ntrl - vars)/2) / 2;
  dterm_unc    = (ln2 - log(Ntrl-1)) / 2;
  bias_unc     = Xdim'*dterm_unc + sum(psiterms_unc);
  
  dterm_cond    = (ln2 - log(Ntrl_y-1)) / 2;
  psiterms_cond = zeros(1,Ym);
  for vi=vars
    idx = (Ntrl_y-vi);
    psiterms_cond = psiterms_cond + psi(idx/2);
  end
  bias_cond = Xdim*dterm_cond + (psiterms_cond/2);

  Hunc  = Hunc  - bias_unc;
  Hcond = Hcond - ones(Nvec,1)*bias_cond;
end

% class weights
w = Ntrl_y ./ Ntrl;

% compute mutual information
I = Hunc - Hcond*w';

% convert to bits
I = I / ln2;

function [Xcen, class_means] = removeclassmeans(X, design)
[Ntrl, Nvec, Ndim] = size(X);
Xcen = X(:,:);
class_means = zeros(Nvec*Ndim,size(design,2));
for k = 1:size(design,2)
  sel = design(:,k);
  tmp = Xcen(sel,:);
  class_means(:,k) = mean(tmp,1);
  Xcen(sel,:) = bsxfun(@minus,tmp,class_means(:,k).');
end
Xcen = reshape(Xcen,[Ntrl Nvec Ndim]);

function Y = indexed2boolean(X)
uX = unique(X);
Y  = false(numel(X),numel(uX));
for k = 1:size(Y,2)
  Y(X==uX(k),k) = true;
end

