function [source] = estimate_fwhm2(source, maxdist)

% ESTIMATE_FWHM2(SOURCE, MAXDIST)
%
% This function computes the Gaussian fwhm of the spatial filters, according to
% least-squares Gaussian fit including data points up until MAXDIST from the
% locations of interest.
% 
% This function can only deal with scalar filters.


if nargin<2, maxdist = 2.5; end % maxdist should be in units of the pos in source
if isempty(maxdist), maxdist = inf; end

if islogical(source.inside)
  inside = find(source.inside);
else
  inside = source.inside;
end
ninside = numel(inside);

if ~isfield(source.avg, 'filter')
  ft_error('the input should contain spatial filters in');
end

nchan   = size(source.avg.filter{inside(1)},2);
ndir    = size(source.avg.filter{inside(1)},1);
if ndir~=1, 
  ft_error('only scalar filters are allowed as input');
end

%get filters and positions
filter = cat(1,source.avg.filter{inside});
pos    = source.pos(inside,:);

%get the filter correlation matrix
Cmat   = filter*filter';
Cmat   = abs(Cmat)./sqrt(diag(Cmat)*diag(Cmat)');

fwhm    = zeros(size(source.pos,1),1);
onesvec = ones(ninside,1);
for k = 1:ninside
 d   = sqrt(sum( (pos-pos(k*onesvec,:)).^2, 2));
 sel = d<=maxdist;
 s   = gaussfit(Cmat(sel,k)',d(sel)');
 fwhm(inside(k)) = s;
end
source.fwhm = fwhm;

% fwhm    = zeros(size(source.pos,1),3,3);
% onesvec = ones(ninside,1);
% for k = 1:ninside
%   dpos = pos-pos(k*onesvec,:);
%   d   = sqrt(sum(dpos.^2, 2));
%   sel = d<=maxdist;
%   [s,~] = gaussfit3D(dpos(sel,:), Cmat(sel,k));
%   fwhm(inside(k),:,:) = s;
% end
% source.fwhm = fwhm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fwhm] = gaussfit(dat, design)

%this function performs a least-squares gaussian fit

%create independent variable
ivar = design.^2;
dat  = log(dat+eps);
beta = dat*ivar'*pinv(ivar*ivar');

sigma = sqrt(-0.5./beta(:,1));
fwhm  = 2.*sqrt(2.*log(2)).*sigma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [invsigma, R] = gaussfit3D(x, y)

% function for a linear least squares fit of a Gaussian function in 3D
% (with a peak of 1 at (0,0,0)).
%
% use as 
%   sigma = gaussfit3D(x, y)
%    
% x = Nx3 and y = Nx1, sigma is the covariance describing the gaussian

% create design matrix based in the independent variable x
design = [x(:,1).^2 x(:,2).^2 x(:,3).^2 2.*(x(:,1).*x(:,2)) 2.*(x(:,1).*x(:,3)) 2.*(x(:,2).*x(:,3))];
%design = design-repmat(mean(design,1),[size(design,1) 1]);
design = cat(2, design, ones(size(design,1),1));

% log-transform the dependent variable y
dat = -2.*log(y+eps);

% regression
beta = design\dat;

% residuals
res  = dat - design*beta;
R    = 1-sum(res.^2)./sum(dat.^2);

% create output
invsigma = [beta(1) beta(4) beta(5);beta(4) beta(2) beta(6);beta(5) beta(6) beta(3)];
%sigma = inv(invsigma);

% by construction the exponent to the gaussian looks like this
%
% exp( -0.5.*(-x'*invsigma*x) )
%
% -x'*invsigma*x = [x1 x2 x3]*[a d e [x1  
%                              d b f  x2
%                              e f c] x3] = 
%
% a*x1^2+b*x2^2+c*x3^2+2d*(x1x2)+2e*(x1x3)+2f*(x2x3)
%
% the variables a through f correspond with the ordered beta weights

