%   CPD_PLOT(X, Y, C); plots 2 data sets. Works only for 2D and 3D data sets.
%
%   Input
%   ------------------ 
%   X           Reference point set matrix NxD;
%   Y           Current postions of GMM centroids;
%   C           (optional) The correspondence vector, such that Y corresponds to X(C,:) 
%
%   See also CPD_REGISTER.

% Copyright (C) 2007 Andriy Myronenko (myron@csee.ogi.edu)
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

function cpd_plot_iter(X, Y, C)

if nargin<2, error('cpd_plot.m error! Not enough input parameters.'); end;
[m, d]=size(Y);

if d>3, error('cpd_plot.m error! Supported dimension for visualizations are only 2D and 3D.'); end;
if d<2, error('cpd_plot.m error! Supported dimension for visualizations are only 2D and 3D.'); end;

% for 2D case
if d==2,
   plot(X(:,1), X(:,2),'r*', Y(:,1), Y(:,2),'bo'); %axis off; axis([-1.5 2 -1.5 2]);
else
% for 3D case
hold on
   plot3(X(:,1),X(:,2),X(:,3),'.','color',[255 223 196]/255)
   plot3(Y(:,1),Y(:,2),Y(:,3),'bo'); % title('X data (red). Y GMM centroids (blue)');set(gca,'CameraPosition',[15 -50 8]);
hold off
end

% plot correspondences
if nargin>2,
    hold on;
    if d==2,
        for i=1:m,
            plot([X(C(i),1) Y(i,1)],[X(C(i),2) Y(i,2)]);
        end
    else
        for i=1:m,
            plot3([X(C(i),1) Y(i,1)],[X(C(i),2) Y(i,2)],[X(C(i),3) Y(i,3)]);
        end
    end
    hold off;
end

drawnow;
