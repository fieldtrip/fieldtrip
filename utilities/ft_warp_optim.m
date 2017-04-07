function [result, M] = ft_warp_optim(input, target, method)

% FT_WARP_OPTIM determine intermediate positions using warping (deformation)
% the input cloud of points is warped to match the target.
% The strategy is to start with simpelest linear warp, followed by a more
% elaborate linear warp, which then is followed by the nonlinear warps up
% to the desired order.
%
% [result, M] = ft_warp_pnt(input, target, method)
%     input          contains the Nx3 measured 3D positions
%     target         contains the Nx3 template 3D positions
%     method         should be any of 
%                     'rigidbody'
%                     'globalrescale'
%                     'traditional' (default)
%                     'nonlin1'
%                     'nonlin2'
%                     'nonlin3'
%                     'nonlin4'
%                     'nonlin5'
%
% The default warping method is a traditional linear warp with individual
% rescaling in each dimension. Optionally you can select a nonlinear warp
% of the 1st (affine) up to the 5th order.
%
% When available, this function will use the MATLAB optimization toolbox.
%
% See also FT_WARP_APPLY, FT_WARP_ERRROR

% Copyright (C) 2000-2013, Robert Oostenveld
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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

% this can be used for printing detailled user feedback
fb = false;

if nargin<3
  method='traditional';
end

pos1 = input;
pos2 = target;

% The ft_warp_error function might be located in the private subdirectory fo
% fieldtrip, i.e. only available to functions in the FieldTrip toolbox.
% The following line ensures that the function can also be found by the
% feval that is executed by the optimalization toolbox.
errorfun = str2func('ft_warp_error');

% determine whether the MATLAB Optimization toolbox is available and can be used
if ft_hastoolbox('optim')
  optimfun = @fminunc;
else
  optimfun = @fminsearch;
end

% set the options for the minimalisation routine
if isequal(optimfun,  @fminunc)
  options  = optimset('fminunc');
  options  = optimset(options, 'Display', 'off');
  options  = optimset(options, 'MaxIter', 1500);
  % options  = optimset(options, 'MaxFunEvals', '1000*numberOfVariables');
  options  = optimset(options, 'TolFun', 1e-4);
  options  = optimset(options, 'LargeScale', 'off');
elseif isequal(optimfun, @fminsearch)
  options  = optimset('fminsearch');
  options  = optimset(options, 'Display', 'off');
  options  = optimset(options, 'MaxIter', 4500);
else
  warning('unknown optimization function "%s", using default parameters', func2str(optimfun));
end

if fb; fprintf('distance = %f\n', errorfun([0 0 0 0 0 0], pos1, pos2, 'rigidbody')); end

if isempty(method)
  error('incorrect warping method specified');
end

% the warp is done in steps, starting simple and progressively getting more complex
level = find(strcmp(method, {
  'rigidbody'         % 1
  'globalrescale'     % 2
  'traditional'       % 3
  'nonlin1'           % 4
  'nonlin2'           % 5
  'nonlin3'           % 6
  'nonlin4'           % 7
  'nonlin5'           % 8
  }));

if isempty(level)
  error('incorrect warping method specified');
end

if level>=1
  % do a rigid-body transformation (6 parameters)
  if fb; disp('rigidbody...'); end
  ri = [0 0 0 0 0 0];
  rf = optimfun(errorfun, ri, options, pos1, pos2, 'rigidbody');
  if fb; fprintf('distance = %f\n', errorfun(rf, pos1, pos2, 'rigidbody')); end
end

if level>=2
  % do a rigid-body + global rescaling transformation (7 parameters)
  if fb; disp('rigidbody + global rescaling...'); end
  gi = [rf 1];
  gf = optimfun(errorfun, gi, options, pos1, pos2, 'globalrescale');
  if fb; fprintf('distance = %f\n', errorfun(gf, pos1, pos2, 'globalrescale')); end
end

if level>=3
  % do a rigid-body + individual rescaling transformation (9 parameters)
  if fb; disp('rigidbody + individual rescaling...'); end
  ti = [gf gf(7) gf(7)];
  tf = optimfun(errorfun, ti, options, pos1, pos2, 'traditional');
  if fb; fprintf('distance = %f\n', errorfun(tf, pos1, pos2, 'traditional')); end
end

if level>=4
  % do a first order nonlinear transformation,
  if fb; disp('1st order nonlinear...'); end
  e1i = traditional(tf);
  e1i = [e1i(1:3,4) e1i(1:3,1:3)];  % reshuffle from homogenous into nonlinear
  e1f = optimfun(errorfun, e1i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', errorfun(e1f, pos1, pos2, 'nonlinear')); end
end

if level>=5
  % do a second order nonlinear transformation,
  if fb; disp('2nd order nonlinear...'); end
  e2i = [e1f zeros(3,6)];
  e2f = optimfun(errorfun, e2i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', errorfun(e2f, pos1, pos2, 'nonlinear')); end
end

if level>=6
  % do a third order nonlinear transformation,
  if fb; disp('3rd order nonlinear...'); end
  e3i = [e2f zeros(3,10)];
  e3f = optimfun(errorfun, e3i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', errorfun(e3f, pos1, pos2, 'nonlinear')); end
end

if level>=7
  % do a fourth order nonlinear transformation,
  if fb; disp('4th order nonlinear...'); end
  e4i = [e3f zeros(3,15)];
  e4f = optimfun(errorfun, e4i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', errorfun(e4f, pos1, pos2, 'nonlinear')); end
end

if level>=8
  % do a fifth order nonlinear transformation,
  if fb; disp('5th order nonlinear...'); end
  e5i = [e4f zeros(3,21)];
  e5f = optimfun(errorfun, e5i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', errorfun(e5f, pos1, pos2, 'nonlinear')); end
end

% return the estimated parameters of the highest level warp
% and compute the warped points
switch level
  case 1
    M = rf;
  case 2
    M = gf;
  case 3
    M = tf;
  case 4
    M = e1f;
  case 5
    M = e2f;
  case 6
    M = e3f;
  case 7
    M = e4f;
  case 8
    M = e5f;
end
result = ft_warp_apply(M, input, method);
