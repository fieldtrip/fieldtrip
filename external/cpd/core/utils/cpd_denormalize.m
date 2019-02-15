%CPD_DENORMALIZE Denormalizes template point set to its original scaling.
%   y=CPD_DENORMALIZE(Y, normal) Denormalizes points in 
%   normalized point set Y to be scaled and shifted
%   back in the reference point set coordinate system.
%
%   Input
%   ------------------ 
%   Y          Normalized point set.
%
%   Output
%   ------------------ 
%   y          denormalized point set back in the reference point set system. 
%
%   normal     structure of scale and shift parameters
%
%   Examples
%   --------
%       x= [1 2; 3 4; 5 6;];
%       y=x;
%       [X, Y0, normal]=cpd_normalize(x,y);
%
%       x2=cpd_denormalize(X, normal);
%       norm(x-x2)
%
%   See also CPD_NORMALIZE, CPD_TRANSFORM, CPD_REGISTER

% Copyright (C) 2006 Andriy Myronenko (myron@csee.ogi.edu)
%
%     This file is part of the Coherent Point Drift (CPD) package.
%
%     The source code is provided under the terms of the GNU General Public License as published by
%     the Free Software Foundation version 2 of the License.
% 
%     CPD package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with CPD package; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

function   Transform =cpd_denormalize(Transform, normal, way);

switch lower(Transform.method)
    case {'rigid','affine'}
        Transform.s=Transform.s*(normal.xscale/normal.yscale);
        Transform.t=normal.xscale*Transform.t+normal.xd'-Transform.s*(Transform.R*normal.yd');
    case 'nonrigid'
        Transform.s=normal.xscale/normal.yscale;
        Transform.t=normal.xd'-Transform.s*normal.yd';
        Transform.W=Transform.W*normal.xscale;   
        Transform.beta=normal.yscale*Transform.beta;       
end







