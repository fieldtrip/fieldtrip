function y = pinvNx2(x)

% PINVNX2 computes a pseudo-inverse of the slices of an Nx2xM real-valued matrix.
% Output has dimensionality 2xNxM. This implementation is generally faster
% than calling pinv in a for-loop, once M > 2 

siz = [size(x) 1];
xtx = zeros([2,2,siz(3:end)]);
xtx(1,1,:,:) = sum(x(:,1,:,:).^2,1);
xtx(2,2,:,:) = sum(x(:,2,:,:).^2,1);
tmp          = sum(x(:,1,:,:).*x(:,2,:,:),1);
xtx(1,2,:,:) = tmp;
xtx(2,1,:,:) = tmp;
ixtx         = inv2x2(xtx);
y            = mtimes2xN(ixtx, permute(x, [2 1 3:ndims(x)]));

function [d] = inv2x2(x)

% INV2X2 computes inverse of matrix x, where x = 2x2xN, using explicit analytic definition

adjx  = [x(2,2,:,:) -x(1,2,:,:); -x(2,1,:,:) x(1,1,:,:)];
denom = x(1,1,:,:).*x(2,2,:,:) - x(1,2,:,:).*x(2,1,:,:);
d     = adjx./denom([1 1],[1 1],:,:);

function [z] = mtimes2xN(x, y)

% MTIMES2XN computes x*y where x = 2x2xM and y = 2xNxM
% and output dimensionatity is 2xNxM 

siz   = size(y);
z     = zeros(siz);

for k = 1:siz(2)
  z(1,k,:,:) = x(1,1,:,:).*y(1,k,:,:) + x(1,2,:,:).*y(2,k,:,:);
  z(2,k,:,:) = x(2,1,:,:).*y(1,k,:,:) + x(2,2,:,:).*y(2,k,:,:);
end