function [warped] = ft_warp_apply(M, input, method, tol)

% FT_WARP_APPLY performs a 3D linear or nonlinear transformation on the input
% coordinates, similar to those in AIR 3.08. You can find technical
% documentation on warping in general at http://bishopw.loni.ucla.edu/AIR3
%
% Use as
%   [warped] = ft_warp_apply(M, input, method, tol)
% where
%   M        vector or matrix with warping parameters
%   input    Nx3 matrix with coordinates
%   warped   Nx3 matrix with coordinates
%   method   string describing the warping method
%   tol      (optional) value determining the numerical precision of the
%             output, to deal with numerical round off imprecisions due to
%             the warping
%
% The methods 'nonlin0', 'nonlin2' ... 'nonlin5' specify a
% polynomial transformation. The size of the transformation matrix
% depends on the order of the warp
%   zeroth order :  1 parameter  per coordinate (translation)
%   first  order :  4 parameters per coordinate (total 12, affine)
%   second order : 10 parameters per coordinate
%   third  order : 20 parameters per coordinate
%   fourth order : 35 parameters per coordinate
%   fifth  order : 56 parameters per coordinate (total 168)
% The size of M should be 3xP, where P is the number of parameters
% per coordinate. Alternatively, you can specify the method to be
% 'nonlinear', where the order will be determined from the size of
% the matrix M.
%
% If the method 'homogeneous' is selected, the input matrix M should be
% a 4x4 homogenous transformation matrix.
%
% If the method 'sn2individual' or 'individual2sn' is selected, the input
% M should be a structure based on nonlinear (warping) normalisation parameters
% created by SPM8 for alignment between an individual structural MRI and the
% template MNI brain.  These options call private functions of the same name.
% M will have subfields like this:
%     Affine: [4x4 double]
%         Tr: [4-D double]
%         VF: [1x1 struct]
%         VG: [1x1 struct]
%      flags: [1x1 struct]
%
% If any other method is selected, it is assumed that it specifies the name of an
% auxiliary function that will, when given the input parameter vector M, return an
% 4x4 homogenous transformation matrix. Supplied functions are 'translate', 'rotate',
% 'scale', 'rigidbody', 'globalrescale', 'traditional', 'affine', 'perspective',
% 'quaternion'.
%
% See also FT_WARP_OPTIM, FT_WARP_ERROR

% Copyright (C) 2000-2016, Robert Oostenveld
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

if nargin<4
  tol = [];
end

if nargin<3 && all(size(M)==4)
  % no specific transformation mode has been selected
  % it looks like a homogenous transformation matrix
  method = 'homogeneous';
elseif nargin<3
  % the default method is 'nonlinear'
  method = 'nonlinear';
end

if size(input,2)==2
  % convert the input points from 2D to 3D representation
  input(:,3) = 0;
  input3d = false;
else
  input3d = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear warping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(method, {'nonlinear', 'nonlin0', 'nonlin1', 'nonlin2', 'nonlin3', 'nonlin4', 'nonlin5'}))
  x = input(:,1);
  y = input(:,2);
  z = input(:,3);
  s = size(M);

  if s(1)~=3
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin0') && s(2)~=1
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin1') && s(2)~=4
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin2') && s(2)~=10
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin3') && s(2)~=20
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin4') && s(2)~=35
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin5') && s(2)~=56
    error('invalid size of nonlinear transformation matrix');
  end

  if s(2)==1
    % this is a translation, which in a strict sense is not the 0th order nonlinear transformation
    xx = M(1,1) + x;
    yy = M(2,1) + y;
    zz = M(3,1) + z;
  elseif s(2)==4
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z;
  elseif s(2)==10
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z;
  elseif s(2)==20
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z;
  elseif s(2)==35
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z + M(1,21)*x.*x.*x.*x + M(1,22)*x.*x.*x.*y + M(1,23)*x.*x.*x.*z + M(1,24)*x.*x.*y.*y + M(1,25)*x.*x.*y.*z + M(1,26)*x.*x.*z.*z + M(1,27)*x.*y.*y.*y + M(1,28)*x.*y.*y.*z + M(1,29)*x.*y.*z.*z + M(1,30)*x.*z.*z.*z + M(1,31)*y.*y.*y.*y + M(1,32)*y.*y.*y.*z + M(1,33)*y.*y.*z.*z + M(1,34)*y.*z.*z.*z + M(1,35)*z.*z.*z.*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z + M(2,21)*x.*x.*x.*x + M(2,22)*x.*x.*x.*y + M(2,23)*x.*x.*x.*z + M(2,24)*x.*x.*y.*y + M(2,25)*x.*x.*y.*z + M(2,26)*x.*x.*z.*z + M(2,27)*x.*y.*y.*y + M(2,28)*x.*y.*y.*z + M(2,29)*x.*y.*z.*z + M(2,30)*x.*z.*z.*z + M(2,31)*y.*y.*y.*y + M(2,32)*y.*y.*y.*z + M(2,33)*y.*y.*z.*z + M(2,34)*y.*z.*z.*z + M(2,35)*z.*z.*z.*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z + M(3,21)*x.*x.*x.*x + M(3,22)*x.*x.*x.*y + M(3,23)*x.*x.*x.*z + M(3,24)*x.*x.*y.*y + M(3,25)*x.*x.*y.*z + M(3,26)*x.*x.*z.*z + M(3,27)*x.*y.*y.*y + M(3,28)*x.*y.*y.*z + M(3,29)*x.*y.*z.*z + M(3,30)*x.*z.*z.*z + M(3,31)*y.*y.*y.*y + M(3,32)*y.*y.*y.*z + M(3,33)*y.*y.*z.*z + M(3,34)*y.*z.*z.*z + M(3,35)*z.*z.*z.*z;
  elseif s(2)==56
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z + M(1,21)*x.*x.*x.*x + M(1,22)*x.*x.*x.*y + M(1,23)*x.*x.*x.*z + M(1,24)*x.*x.*y.*y + M(1,25)*x.*x.*y.*z + M(1,26)*x.*x.*z.*z + M(1,27)*x.*y.*y.*y + M(1,28)*x.*y.*y.*z + M(1,29)*x.*y.*z.*z + M(1,30)*x.*z.*z.*z + M(1,31)*y.*y.*y.*y + M(1,32)*y.*y.*y.*z + M(1,33)*y.*y.*z.*z + M(1,34)*y.*z.*z.*z + M(1,35)*z.*z.*z.*z + M(1,36)*x.*x.*x.*x.*x + M(1,37)*x.*x.*x.*x.*y + M(1,38)*x.*x.*x.*x.*z + M(1,39)*x.*x.*x.*y.*y + M(1,40)*x.*x.*x.*y.*z + M(1,41)*x.*x.*x.*z.*z + M(1,42)*x.*x.*y.*y.*y + M(1,43)*x.*x.*y.*y.*z + M(1,44)*x.*x.*y.*z.*z + M(1,45)*x.*x.*z.*z.*z + M(1,46)*x.*y.*y.*y.*y + M(1,47)*x.*y.*y.*y.*z + M(1,48)*x.*y.*y.*z.*z + M(1,49)*x.*y.*z.*z.*z + M(1,50)*x.*z.*z.*z.*z + M(1,51)*y.*y.*y.*y.*y + M(1,52)*y.*y.*y.*y.*z + M(1,53)*y.*y.*y.*z.*z + M(1,54)*y.*y.*z.*z.*z + M(1,55)*y.*z.*z.*z.*z + M(1,56)*z.*z.*z.*z.*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z + M(2,21)*x.*x.*x.*x + M(2,22)*x.*x.*x.*y + M(2,23)*x.*x.*x.*z + M(2,24)*x.*x.*y.*y + M(2,25)*x.*x.*y.*z + M(2,26)*x.*x.*z.*z + M(2,27)*x.*y.*y.*y + M(2,28)*x.*y.*y.*z + M(2,29)*x.*y.*z.*z + M(2,30)*x.*z.*z.*z + M(2,31)*y.*y.*y.*y + M(2,32)*y.*y.*y.*z + M(2,33)*y.*y.*z.*z + M(2,34)*y.*z.*z.*z + M(2,35)*z.*z.*z.*z + M(2,36)*x.*x.*x.*x.*x + M(2,37)*x.*x.*x.*x.*y + M(2,38)*x.*x.*x.*x.*z + M(2,39)*x.*x.*x.*y.*y + M(2,40)*x.*x.*x.*y.*z + M(2,41)*x.*x.*x.*z.*z + M(2,42)*x.*x.*y.*y.*y + M(2,43)*x.*x.*y.*y.*z + M(2,44)*x.*x.*y.*z.*z + M(2,45)*x.*x.*z.*z.*z + M(2,46)*x.*y.*y.*y.*y + M(2,47)*x.*y.*y.*y.*z + M(2,48)*x.*y.*y.*z.*z + M(2,49)*x.*y.*z.*z.*z + M(2,50)*x.*z.*z.*z.*z + M(2,51)*y.*y.*y.*y.*y + M(2,52)*y.*y.*y.*y.*z + M(2,53)*y.*y.*y.*z.*z + M(2,54)*y.*y.*z.*z.*z + M(2,55)*y.*z.*z.*z.*z + M(2,56)*z.*z.*z.*z.*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z + M(3,21)*x.*x.*x.*x + M(3,22)*x.*x.*x.*y + M(3,23)*x.*x.*x.*z + M(3,24)*x.*x.*y.*y + M(3,25)*x.*x.*y.*z + M(3,26)*x.*x.*z.*z + M(3,27)*x.*y.*y.*y + M(3,28)*x.*y.*y.*z + M(3,29)*x.*y.*z.*z + M(3,30)*x.*z.*z.*z + M(3,31)*y.*y.*y.*y + M(3,32)*y.*y.*y.*z + M(3,33)*y.*y.*z.*z + M(3,34)*y.*z.*z.*z + M(3,35)*z.*z.*z.*z + M(3,36)*x.*x.*x.*x.*x + M(3,37)*x.*x.*x.*x.*y + M(3,38)*x.*x.*x.*x.*z + M(3,39)*x.*x.*x.*y.*y + M(3,40)*x.*x.*x.*y.*z + M(3,41)*x.*x.*x.*z.*z + M(3,42)*x.*x.*y.*y.*y + M(3,43)*x.*x.*y.*y.*z + M(3,44)*x.*x.*y.*z.*z + M(3,45)*x.*x.*z.*z.*z + M(3,46)*x.*y.*y.*y.*y + M(3,47)*x.*y.*y.*y.*z + M(3,48)*x.*y.*y.*z.*z + M(3,49)*x.*y.*z.*z.*z + M(3,50)*x.*z.*z.*z.*z + M(3,51)*y.*y.*y.*y.*y + M(3,52)*y.*y.*y.*y.*z + M(3,53)*y.*y.*y.*z.*z + M(3,54)*y.*y.*z.*z.*z + M(3,55)*y.*z.*z.*z.*z + M(3,56)*z.*z.*z.*z.*z;
  else
    error('invalid size of nonlinear transformation matrix');
  end

  warped = [xx yy zz];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % linear warping using homogenous coordinate transformation matrix
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(method, 'homogenous') || strcmp(method, 'homogeneous')
  if all(size(M)==3)
    % convert the 3x3 homogenous transformation matrix (corresponding with 2D)
    % into a 4x4 homogenous transformation matrix (corresponding with 3D)
    M = [
      M(1,1) M(1,2)  0  M(1,3)
      M(2,1) M(2,2)  0  M(2,3)
      0      0       0  0
      M(3,1) M(3,2)  0  M(3,3)
      ];
  end

  %warped = M * [input'; ones(1, size(input, 1))];
  %warped = warped(1:3,:)';

  % below achieves the same as lines 154-155
  warped = [input ones(size(input, 1),1)]*M(1:3,:)';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % using external function that returns a homogeneous transformation matrix
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif exist(method, 'file') && ~isa(M, 'struct')
  % get the homogenous transformation matrix
  H = feval(method, M);
  warped = ft_warp_apply(H, input, 'homogeneous');

elseif strcmp(method, 'sn2individual') && isa(M, 'struct')
  % use SPM structure with parameters for an inverse warp
  % from normalized space to individual, can be non-linear
  warped = sn2individual(M, input);

elseif strcmp(method, 'individual2sn') && isa(M, 'struct')
  % use SPM structure with parameters for a warp from
  % individual space to normalized space, can be non-linear
  %error('individual2sn is not yet implemented');
  warped = individual2sn(M, input);
else
  error('unrecognized transformation method');
end

if ~input3d
  % convert from 3D back to 2D representation
  warped = warped(:,1:2);
end

if ~isempty(tol)
  if tol>0
    warped = fix(warped./tol)*tol;
  end
end
