function I = mi_mixture_gd(x, y, Ym)
% MI_GD Mutual information (MI) between a Gaussian and a discrete
%         variable in bits calculated from Gaussian mixture
%   I = mi_gd(x,y,Ym) returns the MI between the (possibly multidimensional)
%   Gaussian variable x and the discrete variable y.
%   Rows of x correspond to samples, columns to dimensions/variables. 
%   (Samples first axis)
%   y should contain integer values in the range [0 Ym-1] (inclusive).
%
%   See also: mi_model_gd

% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if ndims(x)~=2
    error('mi_mixture_gd: input arrays should be 2d')
end
if isvector(y)
    y = y(:);
else
    error('mi_mixture_gd: only univariate discrete variable supported');
end

Ntrl = size(x,1);
Nvar = size(x,2);

if size(y,1) ~= Ntrl
    error('mi_mixture_gd: number of trials do not match');
end


Hcond = zeros(1,Ym);
NtrlY = zeros(1,Ym);
m = zeros(Nvar,Ym);
w = zeros(1,Ym);
C = zeros(Nvar,Nvar,Ym);
chC = zeros(Nvar,Nvar,Ym);
for yi=1:Ym
    % class conditional data
    idx = y==(yi-1);
    dc = x(idx,:);
    % class mean
    m(:,yi) = mean(dc);
    % class weight
    NtrlY(yi) = sum(idx);
    w(yi) = NtrlY(yi) / Ntrl;
    
    % copnorm class data
%     cdc = copnorm(dc);
%     cdc = bsxfun(@times, cdc, std(dc));

    cdc = dc;
    cdc = bsxfun(@minus, cdc, m(:,yi)');
    % covariance
    C(:,:,yi) = (cdc'*cdc) / (NtrlY(yi) - 1);
    chC(:,:,yi) = chol(C(:,:,yi));
    % entropy in nats
    Hcond(yi) = sum(log(diag(chC(:,:,yi)))) + 0.5*Nvar*(log(2*pi)+1);
end

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

D = Nvar;
Ds = sqrt(Nvar);
Hmix = 0;
for yi=1:Ym
    Ps = Ds * chC(:,:,yi)';
    % unscented points for this class
    usc = [bsxfun(@plus,Ps,m(:,yi)) bsxfun(@minus,m(:,yi),Ps)];
    
    % class log-likelihoods at unscented points
    log_lik = zeros(Ym,2*Nvar);
    for mi=1:Ym
        % demean points
        dx = bsxfun(@minus,usc,m(:,mi));
        % gaussian likelihood
        log_lik(mi,:) = norm_innerv(dx, chC(:,:,mi)) - Hcond(mi) + 0.5*Nvar;
    end
    % log mixture likelihood for these unscented points
    logmixlik = maxstar(log_lik, w);
    % add to entropy estimate
    Hmix = Hmix + w(yi)*sum(logmixlik);
end
Hmix = -Hmix/(2*D);

% no bias correction
% can correct Hcond, but not Hmix
% good chance bulk of covariance effects cancel
I = (Hmix - w*Hcond') / log(2);



function w = norm_innerv(x, chC)
% normalised innervations
m = (chC')\x;
w = -0.5 *sum(m.*m,1);
