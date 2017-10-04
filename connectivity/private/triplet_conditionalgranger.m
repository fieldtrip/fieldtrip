function Cij = triplet_conditionalgranger(H3,Z3,cmbindx3,H2,Z2,cmbindx2,cmbindx)

% TRIPLET_CONDITIONALGRANGER
% 
% Inputs:
%   H3,Z3: transfer matrix, noise covariance for
%     triplets, 3x3(xtriplet)xnfreq
%   H2,Z2: transfer matrix, noise covariance for
%     duplets,  2x2(xnduplet)xnfreq
%   cmbindx: Nx3 indices determining the output, abc = b->a/c

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

% the ordering for the triplet indices w.r.t. the duplet indices is crucial
% for the interpretation:
%
% abc - ab gives b->a, conditioned on c. If you want a->b, conditioned on
% c, one should also do a bac - ba computation
%
% Here, the cmbindx are taken, and the appropriately ordered corresponding
% cmbindx2 are generated, if necessary, by swapping the order

if nargin<7 || isempty(cmbindx)
    cmbindx = cmbindx3;
end

% FIXME ensure that the number of frequency bins matches

nfreq = size(H3,4);
Cij = zeros(size(cmbindx,1), nfreq);
for k = 1:size(cmbindx,1)
  % triplet
  sel3     = find(all(ismember(cmbindx3, cmbindx(k,:)),2));
  tmp      = cmbindx3(sel3,:);
  reorder3 = [find(tmp==cmbindx(k,1)) find(tmp==cmbindx(k,2)) find(tmp==cmbindx(k,3))];
  
  h123 = reshape(H3(reorder3, reorder3, sel3, :), [3 3 nfreq]);
  z123 = Z3(reorder3, reorder3, sel3);
  
  % normalization matrix 2->1/3
  tmp1 = -z123(2,1)/z123(1,1);
  tmp2 = -z123(3,1)/z123(1,1);
  tmp3 = -(z123(3,2)+tmp2*z123(1,2))/(z123(2,2)+tmp1*z123(1,2));

  p1  = [1    0 0;
         tmp1 1 0;
         tmp2 0 1];
  p2  = [1 0    0;
         0 1    0;
         0 tmp3 1];
     
  P    = p2*p1;
  invP = inv(P);
  
  %bivariate system  
  sel2     = find(all(ismember(cmbindx2, cmbindx(k,[1 3])),2));
  tmp      = cmbindx2(sel2,:);
  reorder2 = [find(tmp==cmbindx(k,1)) find(tmp==cmbindx(k,3))];
    
  h2 = reshape(H2(reorder2, reorder2, sel2, :), [2 2 nfreq]);
  z2 = Z2(reorder2, reorder2, sel2);
  
  Q  = [1                 0;
        -z2(1,2)/z2(1,1) 1];
  numer = z2(1,1);
  
  invh2   = inv2x2(h2);
  
  HH   = mtimes3x3(h123, invP(:,:,ones(size(h123,3),1)));
  B    = mtimes2x2(Q(:,:,ones(size(invh2,3),1)), invh2);
  FF   = shiftdim(B(1,1,:).*HH(1,1,:)+B(1,2,:).*HH(3,1,:),1);
  denom = abs(FF.*conj(FF)).*z123(1,1);
  Cij(k,:) = log(numer./denom);
  
%   % it seems only FF(1,1) is needed (as per the commented for-loop
%   for kk = 1:size(h2,3)
%     HH = h123(:,:,kk)/P;
%     B  = Q/h2(:,:,kk);
%     FF = B(1,:)*HH([1 3],1);
%     denom = abs(FF(1,1)*z123(1,1)*conj(FF(1,1)));
% 
%     Cij(k,kk) = log(numer./denom);
%     
%   end
%   
  
%   for kk = 1:size(h2,3)
%     HH = h123(:,:,kk)/P;
%     Q  = [1                 0;
%           -z2(1,2)/z2(1,1) 1];
%     B  = Q/h2(:,:,kk);
%     BB = [B(1,1) 0 B(1,2);
%           0      1 0;
%           B(2,1) 0 B(2,2)];
%     FF = BB*HH;
%     numer = z2(1,1);
%     denom = abs(FF(1,1)*z123(1,1)*conj(FF(1,1)));
% 
%     Cij(k,kk) = log(numer./denom);
%   end
end
