%CPD_G Construct Gaussian affinity matrix
%   G=CPD_G(x,y,beta) returns Gaussian affinity matrix between x and y data
%   sets. If x=y returns Gaussian Gramm matrix.
%
%   Input
%   ------------------ 
%   x, y       real, full 2-D matrices. Rows represent samples. Columns
%               represent features.
%   
%   beta      std of the G.
%
%   Output
%   ------------------ 
%   G           Gaussian affinity matrix 
%
%   Examples
%   --------
%       x= [1 2; 3 4; 5 6;];
%       beta=2;
%       G=cpd_G(x,x,beta);
%
%   See also CPD_REGISTER.

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

function G=cpd_G(x,y,beta)

if nargin<3, error('cpd_G.m error! Not enough input parameters.'); end;

k=-2*beta^2;
[n, d]=size(x); [m, d]=size(y);

G=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
G=squeeze(sum(G.^2,2));
G=G/k;
G=exp(G);
