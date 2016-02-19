function [h] = procrustes_trans(input,target)

% PROCRUSTES_TRANS returns the homogenous coordinate transformation matrix
% that warps the specified input points to the target points. 
%
% Use as
%   [h] = procrustes_trans(input, target) 
% where
%   input   Nx3 matrix with coordinates
%   target  Nx3 matrix with coordinates
% 
% The algorithm used for the calculation of the rotation matrix is knonwn
% as the Procrustes method. Its use for MEG coordinate transformation has 
% been suggested in Fuchs et al. TBME vol. 42, 1995, p. 416ff.
% 
% See also WARP_OPTIM, HEADCOORDINATES

% Copyright (C) 2010, Tilmann Sander-Thoemmes
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

ninp = size(input,1);
ntarg = size(target,1);

% do basic checks
if ninp ~= ntarg,     
    error('you must specify same number of points for input and target');
end
if  ninp < 3,     
    error('you must specify at least three points for matching');
end

% calculate the center fo gravity
ctr_grav_inp = mean(input);
ctr_grav_targ = mean(target);

% subtract center of gravity from points: set of directions is obtained
ipnt=input-repmat(ctr_grav_inp,ninp,1);
tpnt=target-repmat(ctr_grav_targ,ninp,1);

% Procrustes: Calculate the best rotation matrix that transforms
% input to target directions in the least squares sense
P = ipnt'*tpnt;       
[U,S,V] = svd(P);
% U*S*V'

rotm = V*U';
% det(rotm)
% Correcting for negative determinant: See Wikipedia on Orthogonal 
% Procrustes Problem
if det(rotm)<0
    mirror=eye(3);
    mirror(3,3) = -1;
    rotm = V*mirror*U';
end

% With the rotation matrix the translation can be calculated
trans = ctr_grav_targ'-rotm*ctr_grav_inp';

% compute the full homogenous transformation matrix 
h = eye(4);
h(1:3,1:3) = rotm;
h(1:3,4) = trans';  

