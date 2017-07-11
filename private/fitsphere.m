function [C,R] = fitsphere(pnt)

% FITSPHERE fits the centre and radius of a sphere to a set of points
% using Taubin's method.
%
% Use as
%       [center,radius] = fitsphere(pnt)
% where
%   pnt     = Nx3 matrix with the Carthesian coordinates of the surface points
% and
%   center  = the center of the fitted sphere
%   radius  = the radius of the fitted sphere

% Copyright (C) 2009, Jean Daunizeau (for SPM)
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

x = pnt(:,1);
y = pnt(:,2);
z = pnt(:,3);

% Make sugary one and zero vectors
l = ones(length(x),1);
O = zeros(length(x),1);

% Make design mx
D = [(x.*x + y.*y + z.*z) x y z l];

Dx = [2*x l O O O];
Dy = [2*y O l O O];
Dz = [2*z O O l O];

% Create scatter matrices
M = D'*D;
N = Dx'*Dx + Dy'*Dy + Dz'*Dz;

% Extract eigensystem
[v, evalues] = eig(M);
evalues = diag(evalues);
Mrank = sum(evalues > eps*5*norm(M));

if (Mrank == 5)
  % Full rank -- min ev corresponds to solution
  % Minverse = v'*diag(1./evalues)*v;
  [v,evalues] = eig(inv(M)*N);
  [dmin,dminindex] = max(diag(evalues));
  pvec = v(:,dminindex(1))';
else
  % Rank deficient -- just extract nullspace of M
  % pvec = null(M)';  % this does not work reliably because of inconsistent rank definition
  pvec = v(:,evalues <= eps*5*norm(M))';
  [m,n] = size(pvec);
  if m > 1
    pvec = pvec(1,:);
  end
end

if isempty(pvec)
  ft_warning('was not able to fit a sphere to the surface points');
   C = [NaN NaN NaN];
   R = Inf;
else
   % Convert to (R,C)
   C = -0.5*pvec(2:4) / pvec(1);
   R = sqrt(sum(C*C') - pvec(5)/pvec(1));
end


% if nargout == 1,
%   if pvec(1) < 0
%     pvec = -pvec;
%   end
%   C = pvec;
% else
%   C = -0.5*pvec(2:4) / pvec(1);
%   R = sqrt(sum(C*C') - pvec(5)/pvec(1));
% end
