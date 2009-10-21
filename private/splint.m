function [V2, L2, L1] = splint(elc1, V1, elc2)

% SPLINT computes the spherical spline interpolation and the surface laplacian
% of an EEG potential distribution
% 
% Use as
%   [V2, L2, L1] = splint(elc1, V1, elc2)
% where
%   elc1	electrode positions where potential is known
%   elc2	electrode positions where potential is not known
%   V1		known potential
% and
%   V2		potential at electrode locations in elc2
%   L2		laplacian of potential at electrode locations in elc2
%   L1		laplacian of potential at electrode locations in elc1
%
% See also LAPINT, LAPINTMAT, LAPCAL

% This implements
%   F. Perrin, J. Pernier, O. Bertrand, and J. F. Echallier.
%   Spherical splines for scalp potential and curernt density mapping.
%   Electroencephalogr Clin Neurophysiol, 72:184-187, 1989.
% including their corrections in 
%   F. Perrin, J. Pernier, O. Bertrand, and J. F. Echallier.
%   Corrigenda: EEG 02274, Electroencephalography and Clinical
%   Neurophysiology 76:565.

% Copyright (C) 2003, Robert Oostenveld
% 
% $Log: splint.m,v $
% Revision 1.7  2006/09/07 12:32:53  roboos
% minor change in documentation
%
% Revision 1.6  2004/06/28 07:37:40  roberto
% tested and improved computation of legendre polynomials
%
% Revision 1.5  2004/01/27 12:43:30  roberto
% re-included splint_gh as subfunction and made it adaptive to the number of electrodes
% added laplacian L1 on elc1 to outpu of function
%
% Revision 1.4  2003/06/03 08:29:39  roberto
% *** empty log message ***
%
% Revision 1.3  2003/04/10 09:48:13  roberto
% fixed bug in construction of matrix H (removes the need for avgref)
% fixed minor bug in gh() subfunction (preallocation of space)
% changed from SVD to Gaussian Elimination to solve for coefficients
%
% Revision 1.2  2003/03/11 14:45:37  roberto
% updated help and copyrights
%

N = size(elc1,1);	% number of known electrodes
M = size(elc2,1);	% number of unknown electrodes
T = size(V1,2);		% number of timepoints in the potential
Z = V1;			% potential on known electrodes, can be matrix

% remember the actual size of the sphere
sphere1_scale = mean(sqrt(sum(elc1.^2,2)));
sphere2_scale = mean(sqrt(sum(elc2.^2,2)));

% scale all electrodes towards a unit sphere
elc1 = elc1 ./ repmat(sqrt(sum(elc1.^2,2)), 1, 3);
elc2 = elc2 ./ repmat(sqrt(sum(elc2.^2,2)), 1, 3);

% compute cosine of angle between all known electrodes
CosEii = zeros(N, N);
for i=1:N
  for j=1:N
    CosEii(i,j) = 1 - sum((elc1(i,:) - elc1(j,:)).^2)/2;
  end
end

% compute matrix G of Perrin 1989
[gx, hx] = gh(CosEii);
G = gx;

% Note that the two simultaneous equations (2) of Perrin
%   G*C + T*c0 = Z
%   T'*C = 0
% with C = [c1 c2 c3 ...] can be rewritten into a single linear system
%   H * [c0 c1 c2 ... cN]' = [0 z1 z2 ... zN]'
% with
%   H = [ 0   1   1   1   1   1   1   .  ]
%       [ 1  g11 g12  .   .   .   .   .  ]
%       [ 1  g21 g22  .   .   .   .   .  ]
%       [ 1  g31 g32  .   .   .   .   .  ]
%       [ 1  g41 g42  .   .   .   .   .  ]
%       [ .   .   .   .   .   .   .   .  ]
% that subsequently can be solved using a singular value decomposition
% or another linear solver

% rewrite the linear system according to the comment above and solve it for C
% different from Perrin, this solution for C includes the c0 term
H = ones(N+1,N+1); H(1,1) = 0; H(2:end,2:end) = G;
% C = pinv(H)*[zeros(1,T); Z];         % solve with SVD
C = H \ [zeros(1,T); Z];               % solve by Gaussian elimination

% compute surface laplacian on all known electrodes by matrix multiplication
L1 = hx * C(2:end, :);

% undo the initial scaling of the sphere to get back to real units for the laplacian
L1 = L1 / (sphere1_scale^2);

if all(size(elc1)==size(elc2)) & all(elc1(:)==elc2(:))
  warning('using shortcut for splint');
  % do not recompute gx and hx
else
  % compute cosine of angle between all known and unknown electrodes
  CosEji  = zeros(M, N);
  for j=1:M
    for i=1:N
      CosEji(j,i) = 1 - sum((elc2(j,:) - elc1(i,:)).^2)/2;
    end
  end
  [gx, hx] = gh(CosEji);
end
% compute interpolated potential on all unknown electrodes by matrix multiplication
V2 = [ones(M,1) gx] * C;

% compute surface laplacian on all unknown electrodes by matrix multiplication
L2 = hx * C(2:end, :);

% undo the initial scaling of the sphere to get back to real units for the laplacian
L2 = L2 / (sphere2_scale^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this subfunction implements equations (3) and (5b) of Perrin 1989 
% for simultaneous computation of multiple values
% also implemented as MEX file, but without the adaptive N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gx, hx] = gh(x)
M = 4;          % constant in denominator
% N is the number of terms for series expansion
% N=9 works fine for 32 electrodes, but for 128 electrodes it should be larger
if max(size(x))<=32
  N = 9;
elseif max(size(x))<=64
  N = 14;
elseif max(size(x))<=128
  N = 20;
else
  N = 32;
end
p  = zeros(1,N);
gx = zeros(size(x));
hx = zeros(size(x));
x(find(x>1)) = 1;       % to avoid rounding off errors
x(find(x<-1)) = -1;     % to avoid rounding off errors
% using Matlab function to compute legendre polynomials
% P = zeros(size(x,1), size(x,2), N);
% for k=1:N
%  tmp = legendre(k,x);
%  P(:,:,k) = squeeze(tmp(1,:,:));
% end
if (size(x,1)==size(x,2)) & all(all(x==x'))
  issymmetric = 1;
else
  issymmetric = 0;
end
for i=1:size(x,1)
  if issymmetric
    jloop = i:size(x,2);
  else
    jloop = 1:size(x,2);
  end
  for j=jloop
    % using faster "Numerical Recipies in C" plgndr function (mex file)
    for k=1:N
      p(k) = plgndr(k,0,x(i,j));
    end
    % p = squeeze(P(i, j ,:))';
    gx(i,j) =  sum((2*(1:N)+1) ./ ((1:N).*((1:N)+1)).^M     .* p) / (4*pi);
    hx(i,j) = -sum((2*(1:N)+1) ./ ((1:N).*((1:N)+1)).^(M-1) .* p) / (4*pi);
    if issymmetric
      gx(j,i) = gx(i,j);
      hx(j,i) = hx(i,j);
    end
    % fprintf('computing laplacian %f%% done\n', 100*(((i-1)*size(x,2)+j) / prod(size(x))));
  end
end
