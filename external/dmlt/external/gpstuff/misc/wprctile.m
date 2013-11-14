function y = wprctile(x, p, w)
% WPRCTILE Percentiles of a weighted sample.
%
%    Description
%      Y = PRCTILE(X, P, W) returns percentiles of the values in X. P
%       is a scalar or a vector of percent values. W is a vector of
%       unnormalized weights for samples. Length of W has to be same
%       as length of X. X need to be a a vector. Y is the same size as
%       P, and Y(i) contains the P(i)-th percentile.
%
%      Example
%       y = prctile(x,50,w); % the median of x given sample weights w
%
%    See also wmean, prctile

% BUGS: Accepts only vector valued X

% Copyright (c) 2000-2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

x=sort(x);
p=p./100;
y=p;
ww=cumsum(w);ww=ww./ww(end);
for j=1:length(p)
  wi=min(find(ww>=p(j)));
  if wi==1
    y(j)=x(1);
  else
    w1=ww(wi-1);x1=x(wi-1);
    y(j)=x1+(x(wi)-x1).*(p(j)-w1)./(ww(wi)-w1);
  end
end
