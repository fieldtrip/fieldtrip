function [y,w] = ccabss(x)
% CCABSS -  Blind Source Separation by Canonical Correlation Analysis
%
% Y = CCABSS(X) is the BSS of X=A*S where S is a set of unknown source signals
% and A is an unknown mixing matrix. The components in S are supposed to
% be independent. Y is an estimate of S appart from permutation and scaling.
% For mixed 1-D signals, X is 2-D. The first index refer to the different
% components and the second index refers to the signal parameter (e.g. time)
% For mixed images, X is 3-D where the first index refers to the different 
% mixed images and the second and third indeces are the spatial coordinates.
%
% [Y W] = CCABSS(X) also gives the 'de-mixing' matrix W, such that Y = W'*X.
%
% © 2000 Magnus Borga


switch ndims(x)
 case 2 % 1D signals
  spatial_mode = 0;
  A = x(:,2:end-1);
  B = conv2(x,[1 0 1],'valid'); % Temporal correlation
  [wa wb r] = cc(A,B); % CCA
  y = wa'*x;
 case 3 % 2D signals
  spatial_mode = 1;
  x_size = size(x);
  im_size = x_size(2)*x_size(3);
  ab_size = (x_size(2)-2)*(x_size(3)-2);
  for k = 1:x_size(1) % Flatten 2D-signals after convolution
    X(k,:) = reshape(x(k,:,:),1,im_size);
    A(k,:) = reshape(x(k,2:end-1,2:end-1),1,ab_size);
    B(k,:) = reshape(conv2(squeeze(x(k,:,:)),[0 1 0;1 0 1;0 1 0],'valid'),1,ab_size);
  end
    [wa wb r] = cc(A,B); % CCA
  for k = 1:x_size(1) % Flatten 2D-signals after convolution
    y(k,:,:) = reshape(wa(:,k)'*X,x_size(2),x_size(3));
  end
 otherwise, error('x must be 2- or 3-dimensional.')
end

if nargout > 1
  w = wa;
end

% ------------
% --- CCA ----
% ------------

function [Wx, Wy, r] = cc(X,Y)

% --- Calculate covariance matrices ---

z = [X;Y];
C = cov(z.');
sx = size(X,1);
sy = size(Y,1);
Cxx = C(1:sx, 1:sx) + 10^(-8)*eye(sx);
Cxy = C(1:sx, sx+1:sx+sy);
Cyx = Cxy';
Cyy = C(sx+1:sx+sy, sx+1:sx+sy) + 10^(-8)*eye(sy);
invCyy = inv(Cyy);

% --- Calcualte Wx and r ---

[Wx,r] = eig(inv(Cxx)*Cxy*invCyy*Cyx); % Basis in X
r = sqrt(real(r));      % Canonical correlations

% --- Sort correlations ---

V = fliplr(Wx);     % reverse order of eigenvectors
r = flipud(diag(r));    % extract eigenvalues anr reverse their orrer
[r,I]= sort((real(r))); % sort reversed eigenvalues in ascending order
r = flipud(r);      % restore sorted eigenvalues into descending order
for j = 1:length(I)
  Wx(:,j) = V(:,I(j));  % sort reversed eigenvectors in ascending order
end
Wx = fliplr(Wx);    % restore sorted eigenvectors into descending order

% --- Calcualte Wy  ---

Wy = invCyy*Cyx*Wx;     % Basis in Y
Wy = Wy./repmat(sqrt(sum(abs(Wy).^2)),sy,1); % Normalize Wy
