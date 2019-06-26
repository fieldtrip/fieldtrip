function [varargout] = meg_leadfield1(varargin)

% MEG_LEADFIELD1 magnetic leadfield for a dipole in a homogenous sphere
%
% [lf] = meg_leadfield1(R, pos, ori)
%
% with input arguments
%   R         position dipole
%   pos     position magnetometers
%   ori     orientation magnetometers
%
% The center of the homogenous sphere is in the origin, the field
% of the dipole is not dependent on the sphere radius.
%
% This function is also implemented as MEX file.

% adapted from Luetkenhoener, Habilschrift '92
% optimized for speed using temporary variables
% the mex implementation is a literary copy of this

% Copyright (C) 2002-2008, Robert Oostenveld
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

% compile the missing mex file on the fly
% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);
  mex(mexsrc);
  cd(pwdir);
  success = true;

catch
  % compilation failed
  disp(lasterr);
  error('could not locate MEX file for %s', mexname);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was just created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE FOLLOWING CODE CORRESPONDS WITH THE ORIGINAL IMPLEMENTATION
% function [lf] = meg_leadfield1(R, Rm, Um);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Nchans = size(Rm, 1);
%
% lf = zeros(Nchans,3);
%
% tmp2 = norm(R);
%
% for i=1:Nchans
%   r = Rm(i,:);
%   u = Um(i,:);
%
%   tmp1 = norm(r);
%   % tmp2 = norm(R);
%   tmp3 = norm(r-R);
%   tmp4 = dot(r,R);
%   tmp5 = dot(r,r-R);
%   tmp6 = dot(R,r-R);
%   tmp7 = (tmp1*tmp2)^2 - tmp4^2;  % cross(r,R)^2
%
%   alpha = 1 / (-tmp3 * (tmp1*tmp3+tmp5));
%   A = 1/tmp3 - 2*alpha*tmp2^2 - 1/tmp1;
%   B = 2*alpha*tmp4;
%   C = -tmp6/(tmp3^3);
%
%   if tmp7<eps
%     beta = 0;
%   else
%     beta = dot(A*r + B*R + C*(r-R), u)/tmp7;
%   end
%
%   lf(i,:) = cross(alpha*u  + beta*r, R);
% end
% lf = 1e-7*lf; % multiply with u0/4pi
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fast cross product
% function [c] = cross(a,b);
% c = [a(2)*b(3)-a(3)*b(2) a(3)*b(1)-a(1)*b(3) a(1)*b(2)-a(2)*b(1)];
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fast dot product
% function [c] = dot(a,b);
% c = sum(a.*b);
%
