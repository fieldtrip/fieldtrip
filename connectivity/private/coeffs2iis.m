function iis = coeffs2iis(A,C)

% COEFFS2IIS computes the instantaneous interaction strength based on the
% MVAR-coefficients and a noise covariance matrix. The underlying
% assumption is that the MVAR-models have been fitted in a bivariate
% fashion. It uses the definition according to Vinck et al. Neuroimage 2015
% (108).
%
% Input data:
%   A = 2x2xncmbxorder, matrix with MVAR-coefficients
%   C = 2x2xncmb      , covariance matrices of the noise

siz  = [size(C) 1];
ncmb = siz(3);

P1 = repmat(eye(2), [1 1 ncmb]);
P2 = repmat(eye(2), [1 1 ncmb]);

P1(2,1,:) = -C(1,2,:)./C(2,2,:);
P2(1,2,:) = -C(2,1,:)./C(1,1,:);

B12 = A;
B21 = A;
for k = 1:size(A,4)
  B12(:,:,:,k) = mtimes2x2(P1,A(:,:,:,k));
  B21(:,:,:,k) = mtimes2x2(P2,A(:,:,:,k));
end

num = abs(P1(2,1,:)) + abs(P2(1,2,:));

m12 = zeros(size(num));
m21 = zeros(size(num));
for k = 1:size(B12,4)
  m12 = max(m12, abs(B12(1,2,:,k)));
  m21 = max(m21, abs(B21(2,1,:,k)));
end

iis = shiftdim(num./(m12+m21));
