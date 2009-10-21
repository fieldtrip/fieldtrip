function [px, py, wx, wy, powx, powy] = cancorr(C,x,y,powflag,realflag,trunc)

% CANCORR computes the canonical correlation between multiple variables
%
% Canonical correlation analysis (CCA) is a way of measuring the linear
% relationship between two multidimensional variables. It finds two bases,
% one for each variable, that are optimal with respect to correlations and,
% at the same time, it finds the corresponding correlations.
%
% Use as
%   [px, py, wx, wy] = cancorr(C,x,y)

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: cancorr.m,v $
% Revision 1.3  2009/10/07 09:56:15  jansch
% numerically stable reimplementation according to Borga 1992
%
% Revision 1.2  2007/07/17 09:55:41  jansch
% added possibility to truncate the values on the diagonal, and to output the
% rotated powers
%
% Revision 1.1  2005/05/09 14:17:36  roboos
% new implementation, used in freqdescriptives
%

if nargin<3,
  error('to compute canonical coherence, you need to specify the indices to the cross-spectral density making up the dependent and independent variable');
elseif nargin==3,
  powflag   = 0;
  realflag  = 0;
  trunc     = 0;
elseif nargin==4,
  realflag  = 0;
  trunc     = 0;
elseif nargin==5,
  trunc     = 0;
end

siz = size(C,1);
ind = find(sum(isfinite(C))>0);

xorig = x;
yorig = y;
x     = intersect(ind, x);
y     = intersect(ind, y);

%use the approach as specified by Borga et al 1992: a unified approach to PCA, PLS, MLR and CCA
%that is: solve the generalized eigenvalue problem eig(inv(B)*A) with
%A = [O Cxy; Cyx O], and B = [Cxx O; O Cyy];
indx    = zeros(length(ind), 1);
indy    = zeros(length(ind), 1);
indx(1:length(x))     = 1;
indy(length(x)+1:end) = 1;

Aorig  = [indx*indy' + indy*indx'].*C([x y],[x y]);
Borig  = [indx*indx' + indy*indy'].*C([x y],[x y]);

if realflag,
  A = real(Aorig);
  B = real(Borig);
else
  A = Aorig;
  B = Borig;
end
%if A and B are rank deficient this could lead to non-finite eigenvalues
%[w, p]     = eig(A,B);
[ua,sa,va]  = svds(C(x,x), rank(C(x,x)));
[ub,sb,vb]  = svds(C(y,y), rank(C(y,y)));
U           = [ua' zeros(size(ua,2),size(ub,1)); zeros(size(ub,2),size(ua,1)) ub'];
[w, p]      = eig(U*A*U', U*B*U');

[srt, ind] = sort(diag(abs(p)), 'descend');
w          = w(:, ind);
p          = p(ind, ind);

%eigenvalues come in pairs, with 180-degree ambiguity
nump = ceil(size(p,1)/2);
%wx   = zeros(length(xorig)); wx(find(ismember(xorig,x)), 1:nump) = w(1:length(x),     1:2:end);
%wy   = zeros(length(yorig)); wy(find(ismember(yorig,y)), 1:nump) = w(length(x)+1:end, 1:2:end);
wx = zeros(length(xorig));
wy = zeros(length(yorig));
px   = abs(p(1:2:end, 1:2:end));
py   = abs(p(1:2:end, 1:2:end));

if powflag,
  powx = wx'*B(x,x)*wx;
  powy = wy'*B(y,y)*wy;
end

if 0,
Cxx = C(x,x);
Cyy = C(y,y);
Cxy = C(x,y);
Cyx = C(y,x);

[wx, px] = eig(inv(Cxx)*Cxy*inv(Cyy)*Cyx);
[wy, py] = eig(inv(Cyy)*Cyx*inv(Cxx)*Cxy);

[srtx,indx] = sort(diag(px), 'descend');
[srty,indy] = sort(diag(py), 'descend');

px = px(indx,indx);
wx = wx(:,   indx);
py = py(indy,indy);
wy = wy(:,   indy);

if powflag,
  powx = wx'*Cxx*wx;
  powy = wy'*Cyy*wy;
end
end

if trunc>0 && trunc<1,
  px = px.*double(px>trunc);
  py = py.*double(py>trunc);
elseif trunc>=1,
  px(trunc+1:end,trunc+1:end) = 0;
  py(trunc+1:end,trunc+1:end) = 0;
end
