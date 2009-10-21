function [result, M] = warp_pnt(input, target, method)

% WARP_PNT determine intermediate positions using warping (deformation)
% the input cloud of points is warped to match the target.
% The strategy is to start with simpelest linear warp, followed by a more
% elaborate linear warp, which then is followed by the nonlinear warps up
% to the desired order.
%
% [result, M] = warp_pnt(input, target, method)
%     input          contains the Nx3 measured 3D positions
%     target         contains the Nx3 template 3D positions
%     method         should be empty or any of 'nonlin1', 'nonlin2' ... 'nonlin5'
%
% The default is a traditional linear warp with rescaling in each
% dimension. Optionally you can select a nonlinear warp of the 1st (affine)
% up to the 5th order.
%
% This function depends on the OPTIM and WARPING toolboxes.

% Copyright (C) 2000-2005, Robert Oostenveld
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

% $Log: warp_optim.m,v $
% Revision 1.8  2009/04/16 12:55:15  roboos
% small change in handling of the optimizer function
%
% Revision 1.7  2006/11/23 11:34:51  roboos
% use optimization toolbox if possible, othewise use the standard fminsearch function
%
% Revision 1.6  2006/09/13 09:55:58  roboos
% fixed bug (typo) in rigidbody
%
% Revision 1.5  2006/04/13 12:55:45  roboos
% added a str2func to solve a problem with feval and private
%
% Revision 1.4  2006/04/13 10:50:34  roboos
% renamed calls to warp3d into warp_apply
%
% Revision 1.3  2006/04/13 10:47:39  roboos
% renamed all calls to warpfun into warp_error
%
% Revision 1.2  2006/04/13 10:38:24  roboos
% fixed a problem due to find/strmatch
%
% Revision 1.1  2005/08/15 08:12:40  roboos
% Renamed warp_pnt into warp_optim for consistency with other functions.
% Also changed the code, the subsequent ordering of the simple to
% more complex warps is handled more clean.
%
% Revision 1.4  2005/03/21 15:43:42  roboos
% fixed bug in output for nonlinear warping
% added support for rigidbody or globalrescale warp
%
% Revision 1.3  2004/05/19 09:57:08  roberto
% added GPL copyright statement, added CVS log item
%

global fb;

if nargin<3
  method='traditional';
end

pos1 = input;
pos2 = target;

% The warp_error function might be located in the private subdirectory fo
% fieldtrip, i.e. only available to functions in the fieldtrip toolbox.
% The following line ensures that the function can also be found by the
% feval that is executed by the optimalization toolbox.
warp_error = str2func('warp_error');


% determine whether the Matlab Optimization toolbox is available and can be used
if hastoolbox('optim')
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

if fb; fprintf('distance = %f\n', warp_error([0 0 0 0 0 0], pos1, pos2, 'rigidbody')); end

% the warp is done in steps, starting simple and progressively getting more complex
level = find(strcmp(method, {
  'rigidbody'          % 1
  'globalrescale'     % 2
  'traditional'       % 3
  'nonlin1'           % 4
  'nonlin2'           % 5
  'nonlin3'           % 6
  'nonlin4'           % 7
  'nonlin5'           % 8
  }));

if isempty(method)
  error('incorrect warping method specified');
end

if level>=1
  % do a rigid-body transformation (6 parameters)
  if fb; disp('rigidbody...'); end
  ri = [0 0 0 0 0 0];
  rf = optimfun(warp_error, ri, options, pos1, pos2, 'rigidbody');
  if fb; fprintf('distance = %f\n', warp_error(rf, pos1, pos2, 'rigidbody')); end
end

if level>=2
  % do a rigid-body + global rescaling transformation (7 parameters)
  if fb; disp('rigidbody + global rescaling...'); end
  gi = [rf 1];
  gf = optimfun(warp_error, gi, options, pos1, pos2, 'globalrescale');
  if fb; fprintf('distance = %f\n', warp_error(gf, pos1, pos2, 'globalrescale')); end
end

if level>=3
  % do a rigid-body + individual rescaling transformation (9 parameters)
  if fb; disp('rigidbody + individual rescaling...'); end
  ti = [gf gf(7) gf(7)];
  tf = optimfun(warp_error, ti, options, pos1, pos2, 'traditional');
  if fb; fprintf('distance = %f\n', warp_error(tf, pos1, pos2, 'traditional')); end
end

if level>=4
  % do a first order nonlinear transformation,
  if fb; disp('1st order nonlinear...'); end
  e1i = traditional(tf);
  e1i = [e1i(1:3,4) e1i(1:3,1:3)];	% reshuffle from homogenous into nonlinear
  e1f = optimfun(warp_error, e1i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', warp_error(e1f, pos1, pos2, 'nonlinear')); end
end

if level>=5
  % do a second order nonlinear transformation,
  if fb; disp('2nd order nonlinear...'); end
  e2i = [e1f zeros(3,6)];
  e2f = optimfun(warp_error, e2i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', warp_error(e2f, pos1, pos2, 'nonlinear')); end
end

if level>=6
  % do a third order nonlinear transformation,
  if fb; disp('3rd order nonlinear...'); end
  e3i = [e2f zeros(3,10)];
  e3f = optimfun(warp_error, e3i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', warp_error(e3f, pos1, pos2, 'nonlinear')); end
end

if level>=7
  % do a fourth order nonlinear transformation,
  if fb; disp('4th order nonlinear...'); end
  e4i = [e3f zeros(3,10)];
  e4f = optimfun(warp_error, e4i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', warp_error(e4f, pos1, pos2, 'nonlinear')); end
end

if level>=8
  % do a fifth order nonlinear transformation,
  if fb; disp('5th order nonlinear...'); end
  e5i = [e4f zeros(3,10)];
  e5f = optimfun(warp_error, e5i, options, pos1, pos2);
  if fb; fprintf('distance = %f\n', warp_error(e5f, pos1, pos2, 'nonlinear')); end
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
result = warp_apply(M, input, method);
