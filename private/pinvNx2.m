function y = pinvNx2(x)

% PINVNX2 computes a pseudo-inverse of the M slices of an MxNx2 real-valued matrix.
% Output has dimensionality Mx2xN. This implementation is generally faster
% than calling pinv in a for-loop, once M > 2 

siz = [size(x) 1];
xtx = zeros([siz(1),2,2]);
xtx(:,1,1) = sum(x(:,:,1).^2,2);
xtx(:,2,2) = sum(x(:,:,2).^2,2);
tmp        = sum(x(:,:,1).*x(:,:,2),2);
xtx(:,1,2) = tmp;
xtx(:,2,1) = tmp;
ixtx       = inv2x2(xtx);
y          = mtimes2xN(ixtx, x);

function [d] = inv2x2(x)

% INV2X2 computes inverse of matrix x, where x = 2x2xN, using explicit analytic definition

%adjx  = [x(:,2,2) -x(:,1,2); -x(:,2,1) x(:,1,1)];
d     = x;
denom = x(:,1,1).*x(:,2,2) - x(:,1,2).*x(:,2,1);
%d     = adjx./denom(:,[1 1],[1 1]);
d(:,1,1) = x(:,2,2)./denom;
d(:,1,2) = -x(:,2,1)./denom;
d(:,2,1) = -x(:,1,2)./denom;
d(:,2,2) = x(:,1,1)./denom;

function [z] = mtimes2xN(x, y)

% MTIMES2XN computes x*y where x = Mx2x2 and y = MxNx2
% and output dimensionality is Mx2xN 

siz   = size(y);
z     = zeros(siz([1 3 2]));

for k = 1:siz(2)
  z(:,1,k) = x(:,1,1).*y(:,k,1) + x(:,1,2).*y(:,k,2);
  z(:,2,k) = x(:,2,1).*y(:,k,1) + x(:,2,2).*y(:,k,2);
end
