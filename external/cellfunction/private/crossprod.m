function [y] = crossprod(x,p,flag)

if nargin<3
  % always augment with constant term
  flag = 1;
end
if nargin<2
  % length of smoothing boxcar
  p = 1;
end
flag = double(flag);

% FIXME only works in case of dim=2
n    = size(x);
indx = tril(ones(n(1)))==1;
N    = 0.5*n(1)*(n(1)+1);
y    = zeros(N+flag, n(2));
for k = 1:n(2)
  tmp    = x(:,k)*x(:,k)';
  y(flag+(1:N),k) = tmp(indx);
end

if flag
  y(1,:) = 1;
end

if p>1
  krn = ones(1,p)./p;
  y   = convn(y, krn, 'valid');
end

% function [y] = crossprod(x, ix)
% 
% %FIXME works only in case of dim=2
% n = size(x);
% y = zeros(0.5*n(1)*(n(1)+1), n(2));
% for k = 1:size(ix,1)
%   y(k,:) = x(ix(k,1),:).*x(ix(k,2),:);
% end

