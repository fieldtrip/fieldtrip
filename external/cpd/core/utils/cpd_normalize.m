%CPD_NORMALIZE Normalizes reference and template point sets to have zero
%mean and unit variance.
%   G=CPD_NORMALIZE(x,y) Normalizes x and y to have zero mean and unit
%   variance.
%
%   Input
%   ------------------ 
%   x, y       real, full 2-D matrices. Rows represent samples. Columns
%              represent features. x - reference point set. y - template
%              point set.
%
%   Output
%   ------------------ 
%   X, Y      normalized reference and template point sets. 
%
%   normal     structure that can be used to rescale and shift the point
%              sets back to its original scaling and position (use cpd_denormalize.m)
%
%   Examples
%   --------
%       x= [1 2; 3 4; 5 6;];
%       y=x;
%       [X, Y, normal]=cpd_normalize(x,y);
%
%       x2=cpd_denormalize(X, normal);
%       norm(x-x2)
%
%   See also CPD_DENORMALIZE, CPD_REGISTER

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

function  [X, Y, normal] =cpd_normalize(x,y)

if nargin<2, error('cpd_normalize error! Not enough input parameters.'); end;

n=size(x,1);
m=size(y,1);

normal.xd=mean(x);
normal.yd=mean(y);

x=x-repmat(normal.xd,n,1);
y=y-repmat(normal.yd,m,1);

normal.xscale=sqrt(sum(sum(x.^2,2))/n);
normal.yscale=sqrt(sum(sum(y.^2,2))/m);

X=x/normal.xscale;
Y=y/normal.yscale;



