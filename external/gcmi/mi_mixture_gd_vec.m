function I = mi_mixture_gd_vec(x, y, Ym)
% MI_MIXTURE_GD_VEC Vectorized MI calculation between multiple Gaussian variables 
%         and a common discrete variable in bits from Gaussian mixture. 
%   I = mi_mixture_gd_vec(x,y,Ym) returns the MI between the (possibly multidimensional)
%   Gaussian variables x and the discrete variable y.
%   size(x) = [Ntrl Nvec Ndim]
%   so each output I(i) = mi_mixture_gd(squeeze(x(:,i,:)), y, Ym);
%
%   y should contain integer values in the range [0 Ym-1] (inclusive).
%
%   biascorrect : true / false option (default true) which specifies whether
%   bias correction should be applied to the esimtated MI.
%   demeaned : false / true option (default false) which specifies whether the
%   input data already has zero mean (true if it has been copula-normalized)
%   See also: MI_MIXTURE_GD, MI_MODEL_GD_VEC

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


for vi1=1:Xdim
  % all voxels for this dimension
  x1 = Xcen(:,:,vi1);
  for yi=1:Ym
    tmp = x1(Yhot(:,yi),:);
    Cm(:,vi1,vi1,yi) = sum(tmp.^2);
  end

  for vi2=(vi1+1):Xdim
    x2  = Xcen(:,:,vi2);
    for yi=1:Ym
      tmp = transpose(sum(x1(Yhot(:,yi),:).*x2(Yhot(:,yi),:)));
      Cm(:,vi1,vi2,yi) = tmp;
      Cm(:,vi2,vi1,yi) = tmp;
    end
  end
end


for yi=1:Ym
  Cm(:,:,:,yi) = Cm(:,:,:,yi) / (Ntrl_y(yi) - 1);
  Cm(:,:,:,yi) = vecchol(Cm(:,:,:,yi));
  for vi=1:Xdim
    Hcond(:,yi) = Hcond(:,yi)+log(Cm(:,vi,vi,yi));
  end
end
% class weights
w = Ntrl_y ./ Ntrl;
% normalise entropy properly
Hcond = Hcond + 0.5*Xdim*(log(2*pi)+1);

% mixture entropy via unscented transform
% See:
% Huber, Bailey, Durrant-Whyte and Hanebeck
% "On entropy approximation for Gaussian mixture random vectors"
% http://dx.doi.org/10.1109/MFI.2008.4648062
%
% Goldberger, Gordon, Greenspan
% "An efficient image similarity measure based on approximations of 
% KL-divergence between two Gaussian mixtures"
% http://dx.doi.org/10.1109/ICCV.2003.1238387

D = Xdim;
Ds = sqrt(Xdim);
Hmix = zeros(Nvec,1);
chC = Cm;

m = reshape(class_means,[Nvec Xdim Ym]);
for yi=1:Ym
    % vector transpose
    Ps = Ds * permute(chC(:,:,:,yi), [1 3 2]);
    % unscented points for this class
    usc = cat(3,bsxfun(@plus,Ps,m(:,:,yi)), bsxfun(@minus,m(:,:,yi),Ps));
    
    % class log-likelihoods at unscented points
    log_lik = zeros(Ym,Nvec,2*Xdim);
    for mi=1:Ym
        % demean points
        dx = bsxfun(@minus,usc,m(:,:,mi));
        % gaussian likelihood
        log_lik(mi,:,:) = bsxfun(@minus,norm_innerv(dx, chC(:,:,:,mi)), Hcond(:,mi)) + 0.5*Xdim;
    end
    % log mixture likelihood for these unscented points
    logmixlik = maxstar(log_lik, w);
    % add to entropy estimate
    Hmix = Hmix + w(yi)*sum(logmixlik,2);
end
Hmix = -Hmix/(2*D);

% compute mutual information
I = Hmix - Hcond*w';

% convert to bits
I = I / log(2);

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

function w = norm_innerv(x, chC)
% normalised innervations
Nvec = size(chC,1);
w = zeros(Nvec,size(x,3));
for vi=1:Nvec
  m = (squeeze(chC(vi,:,:))')\(squeeze(x(vi,:,:)));
  w(vi,:) = -0.5 *sum(m.*m,1);
end
