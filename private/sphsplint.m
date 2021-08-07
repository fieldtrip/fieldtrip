function [WVo, WLo] = sphsplint(elc1, elc2, order, degree, lambda)

% SPHSPLINT computes the spherical spline interpolation and the surface
% laplacian of an EEG potential distribution
% 
% Use as
%   [WVo, WLo] = sphsplint(elc1, elc2)
%   [WVo, WLo] = sphsplint(elc1, elc2, order, degree, lambda)
% where
%   elc1    electrode positions where potential is known
%   elc2    electrode positions where potential is not known
% and
%   WVo     filter for the potential at electrode locations in elc2
%   WLo     filter for the laplacian at electrode locations in elc2
%   order   order of splines
%   degree  degree of Legendre polynomials
%   lambda  regularization parameter
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
% Copyright (C) 2021, Ricardo Bruña
% 
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$


% Sets the defaults.
if (nargin < 2), elc2 = elc1; end
if (nargin < 3 || isempty(order)),  order  = 4;    end
if (nargin < 4 || isempty(degree)), degree = [];   end
if (nargin < 5 || isempty(lambda)), lambda = 1e-5; end

% Gets the size of the input set of electrodes.
siz1 = size(elc1, 1);

% Sets the default degree, if no provided.
if isempty(degree)
  if     siz1 <= 32,  degree = 9;
  elseif siz1 <= 64,  degree = 14;
  elseif siz1 <= 128, degree = 20;
  else,               degree = 32;
  end
end


% Creates a single set with all the electrodes.
elcs = unique(cat(1, elc1, elc2), 'rows');
[dum1, ind1] = ismember(elc1, elcs, 'rows');
[dum2, ind2] = ismember(elc2, elcs, 'rows');

% Fits a sphere to the electrodes and centers them.
[center, radius] = fitsphere(elcs);
elcs = elcs - center;

% % Checks the spherical fitting ->  this is already done in the caller
% function
% err    = mean ( sqrt ( sum ( elcs .^2, 2 ) ) ) / radius - 1;
% if err < 0.01
%   fprintf ( 'Perfect spherical fit (residual: %.1f%%).\n', 100 * err );
% elseif err < 0.1
%   fprintf ( 'Good spherical fit (residual: %.1f%%).\n', 100 * err );
% else
%   ft_warning ( 'Bad spherical fit (residual: %.2f%%). The interpolation will be inaccurate.', 100 * err );
% end

% Projects the electrodes over a unit-sphere.
elcs = elcs ./ sqrt(sum(elcs .^ 2, 2));

% Gets the cosine of the angle (projection) between electrodes (Eq. 4).
CosE = min(max(elcs * elcs(ind1, :)', -1), 1);

% Solves the g(x) and h(x) equations (Eqs. 3 and 5b).
[gx, hx] = perrin_gh(CosE, order, degree);

% Gets the g(x) submatrix for the input set alone.
gxii = gx(ind1, :);

% We can rewrite the simultaneous Eq. 2:
%   G * C + T * c0 = Z;
%   T' * C = 0;
% as:
%   H * [ c0 C ]' = [ 0 Z ]'
% with:
%   H = [ 0 1 ]
%       [ 1 G ]

% Builds the composite matrix H (with Tikhonov regularization for g).
H = ones(siz1 + 1);
H(1, 1)         = 0;
H(2:end, 2:end) = gxii + eye(siz1) * lambda;

% Fits the solution for the input electrodes using linear estimation.
iH   = pinv(H);

% Gets the g(x) and h(x) submatrices for the output set.
gxio = gx(ind2, :);
hxio = hx(ind2, :);

% Computes the filter to interpolate the potential in the output set.
WVo  = cat(2, ones(size(gxio, 1), 1), gxio) * iH (:, 2:end);

% Computes the filter to interpolate the Laplacian in the output set.
WLo  = hxio * iH(2:end, 2:end);

% Scales the Laplacian to the original geometrical dimensions.
WLo  = WLo / (radius ^ 2);

function [gx, hx] = perrin_gh(x, order, degree)

% Calculates the Legendre polynomials up to the requested degree.
lpol = mplgndr(degree, 0, x);
lpol = lpol(:, 2:end);

% Gets the g(x) and h(x) functions (Eqs. 3 and 5b).
gx   =  sum ((2 * (1:degree) + 1) ./ ((1:degree) .* ((1:degree) + 1)) .^ order .* lpol, 2) / (4 * pi);
hx   = -sum ((2 * (1:degree) + 1) ./ ((1:degree) .* ((1:degree) + 1)) .^ (order - 1) .* lpol, 2) / (4 * pi);

% Reshapes the g(x) and h(x) functions as matrices.
gx   = reshape(gx, size(x));
hx   = reshape(hx, size(x));
