function [H] = quaternion(q)

% QUATERNION returns the homogenous coordinate transformation matrix corresponding to
% a coordinate transformation described by 7 quaternion parameters.
%
% Use as
%   [H] = quaternion(Q)
% where
%   Q       [q0, q1, q2, q3, q4, q5, q6] vector with parameters. In case
%           of vector of length=6 it is assumed q0 is missing [the case for
%           quanternions from Neuromag MaxFilter output].
%   H       corresponding homogenous transformation matrix
%
<<<<<<< HEAD
% See Elekta/Neuromag MaxFilter manual version 2.2, section "D2 Coordinate Matching",
% page 77 for more details and
=======
% If the input vector has length 6, it is assumed to represent a unit quaternion without scaling.
%
% See Neuromag/Elekta MaxFilter manual version 2.2, section "D2 Coordinate Matching", page 77 for more details and
>>>>>>> cef9aba30fb05e0319a6461d4a246e3bdf321a4b
% https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Conversion_to_and_from_the_matrix_representation

% Copyright (C) 2016, Robert Oostenveld. [!] Edit: 2017-11-27 by MCV.
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

<<<<<<< HEAD
if numel(q)==7
    q0 = q(0+1);
    q1 = q(1+1);
    q2 = q(2+1);
    q3 = q(3+1);
    q4 = q(4+1);
    q5 = q(5+1);
    q6 = q(6+1);
    
elseif numel(q)==6  %First q not present (for MaxFilter quaternions)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q5 = q(5);
    q6 = q(6); 
    
%     q0 = sqrt(1-sum([q1,q2,q3]).^2);  
    
    q0 = sqrt(1-q1^2-q2^2-q3^2);  
=======
if numel(q)==6
  % this is used a lot by the Neuromag/Elekta software, where the first element is left out and a rigid body transformation wothout scaling is used.
  % see also https://github.com/mne-tools/mne-python/blob/maint/0.15/mne/transforms.py#L1137
  q0 = sqrt(1 - q(1)^2 - q(2)^2 - q(3)^2);
  q = [q0 q];
end
>>>>>>> cef9aba30fb05e0319a6461d4a246e3bdf321a4b

    
else
  ft_error('incorrect input vector');
end



R = [
  q0^2+q1^2-q2^2-q3^2  2*(q1*q2-q0*q3)      2*(q1*q3+q0*q2)
  2*(q1*q2+q0*q3)      q0^2+q2^2-q1^2-q3^2  2*(q2*q3-q0*q1)
  2*(q1*q3-q0*q2)      2*(q2*q3+q0*q1)      q0^2+q3^2-q1^2-q2^2
  ];

T = [q4 q5 q6]';

H = eye(4,4);
H(1:3,1:3) = R;
H(1:3,4)   = T;
