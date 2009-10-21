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
% $Log: fitsphere.m,v $
% Revision 1.5  2009/05/05 10:02:42  vlalit
% Typo fix
%
% Revision 1.4  2009/05/01 07:32:27  roboos
% added another solution for the flat surface fitting: now the norm is being used consistently
%
% Revision 1.3  2009/04/30 16:59:36  vlalit
% Bug fix for the problem of too flat mesh surfaces as suggested by Christophe
%
% Revision 1.2  2009/04/01 12:37:12  roboos
% updated documentation
%
% Revision 1.1  2009/04/01 12:16:14  roboos
% new function, see email from Guillaume Flandin from 23 March
%

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
  warning('was not able to fit a sphere to the surface points');
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
