function [mnv,mxv] = spm_minmax(g)
% Compute a suitable range of intensities for VBM preprocessing stuff
% FORMAT [mnv,mxv] = spm_minmax(g)
% g    - array of data
% mnv  - minimum value
% mxv  - maximum value
%
% A MOG with two Gaussians is fitted to the intensities.  The lower
% Gaussian is assumed to represent background.  The lower value is
% where there is a 50% probability of being above background.  The
% upper value is one that encompases 99.5% of the values.
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


d   = [size(g) 1];
mxv = double(max(g(:)));
mnv = double(min(g(:)));
h   = zeros(256,1);
spm_progress_bar('Init',d(3),'Initial histogram','Planes loaded');
sw = warning('off','all');
for i=1:d(3)
    h = h + spm_hist(uint8(round(g(:,:,i)*(255/mxv))),ones(d(1)*d(2),1));
    spm_progress_bar('Set',i);
end;
warning(sw);
spm_progress_bar('Clear');

% Occasional problems with partially masked data because the first Gaussian
% just fits the big peak at zero.  This will fix that one, but cause problems
% for skull-stripped data.
h(1)   = 0;
h(end) = 0;

% Very crude heuristic to find a suitable lower limit.  The large amount
% of background really messes up mutual information registration.
% Begin by modelling the intensity histogram with two Gaussians.  One
% for background, and the other for tissue.
[mn,v,mg] = fithisto((0:255)',h,1000,[1 128],[32 128].^2,[1 1]);
pr        = distribution(mn,v,mg,(0:255)');

%fg = spm_figure('FindWin','Interactive');
%if ~isempty(fg)
%    figure(fg);
%    plot((0:255)',h/sum(h),'b', (0:255)',pr,'r');
%    drawnow;
%end;

% Find the lowest intensity above the mean of the first Gaussian
% where there is more than 50% probability of not being background
mnd       = find((pr(:,1)./(sum(pr,2)+eps) < 0.5) & (0:255)' >mn(1));
if isempty(mnd) || mnd(1)==1 || mn(1)>mn(2),
    mnd = 1;
else
    mnd = mnd(1)-1;
end
mnv       = mnd*mxv/255;

% Upper limit should get 99.5% of the intensities of the
% non-background region
ch  = cumsum(h(mnd:end))/sum(h(mnd:end));
ch  = find(ch>0.995)+mnd;
mxv = ch(1)/255*mxv;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function [mn,v,mg,ll]=fithisto(x,h,maxit,n,v,mg)
% Fit a mixture of Gaussians to a histogram
h   = h(:);
x   = x(:);
sml = mean(diff(x))/1000;
if nargin==4
    mg = sum(h);
    mn = sum(x.*h)/mg;
    v  = (x - mn); v = sum(v.*v.*h)/mg*ones(1,n);
    mn = (1:n)'/n*(max(x)-min(x))+min(x);
    mg = mg*ones(1,n)/n;
elseif nargin==6
    mn = n;
    n  = length(mn);
else
    error('Incorrect usage');
end;

ll = Inf;
for it=1:maxit
    prb  = distribution(mn,v,mg,x);
    scal = sum(prb,2)+eps;
    oll  = ll;
    ll   = -sum(h.*log(scal));
    if it>2 &&  oll-ll < length(x)/n*1e-9
        break;
    end;
    for j=1:n
        p     = h.*prb(:,j)./scal;
        mg(j) = sum(p);
        mn(j) = sum(x.*p)/mg(j);
        vr    = x-mn(j);
        v(j)  = sum(vr.*vr.*p)/mg(j)+sml;
    end;
    mg = mg + 1e-3;
    mg = mg/sum(mg);
end;
%_______________________________________________________________________
%_______________________________________________________________________
function y=distribution(m,v,g,x)
% Gaussian probability density
if nargin ~= 4
    error('not enough input arguments');
end;
x = x(:);
m = m(:);
v = v(:);
g = g(:);
if ~all(size(m) == size(v) & size(m) == size(g))
    error('incompatible dimensions');
end;

for i=1:size(m,1)
    d      = x-m(i);
    amp    = g(i)/sqrt(2*pi*v(i));
    y(:,i) = amp*exp(-0.5 * (d.*d)/v(i));
end;
return;
