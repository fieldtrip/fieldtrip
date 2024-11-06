function [circumference] = ft_headcircumference(cfg, headshape)

% FT_HEADCIRCUMFERENCE determines the head circumference from a triangulated mesh of
% the scalp in the same way as it would be measured using a measuring tape for
% fitting an EEG cap.
%
% Use as
%   circumference = ft_headcircumference(cfg, mesh)
% where the input mesh corresponds to the output of FT_PREPARE_MESH.
%
% The configuration should contain
%   cfg.fiducial.nas   = 1x3 vector with coordinates
%   cfg.fiducial.ini   = 1x3 vector with coordinates
%   cfg.fiducial.lpa   = 1x3 vector with coordinates
%   cfg.fiducial.rpa   = 1x3 vector with coordinates
%   cfg.feedback       = string, can be 'yes' or 'no' for detailed feedback (default = 'yes')
%
% See also FT_ELECTRODEPLACEMENT, FT_PREPARE_MESH, FT_VOLUMESEGMENT, FT_READ_HEADSHAPE

% Copyright (C) 2024, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar mesh
ft_preamble provenance mesh

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% convert to structure
if isnumeric(headshape) && size(headshape,2)==3
  headshape = struct('pos', headshape);
end

% add triangulation
if ~isfield(headshape, 'tri') && isfield(headshape, 'pos')
  prj = elproj(headshape.pos);
  headshape.tri = delaunay(prj(:,1), prj(:,2));
end

% check if the input data is valid for this function
headshape = ft_checkdata(headshape, 'datatype', 'mesh', 'hasunit', 'yes');

% set the defaults
cfg.feedback = ft_getopt(cfg, 'feedback', 'yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is largely shared with the 1020 method in FT_ELECTRODEPLACEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the placement procedure fails if the fiducials coincide with vertices
dist = @(x, y) sqrt(sum(bsxfun(@minus, x, y).^2,2));
tolerance = 0.1 * ft_scalingfactor('mm', headshape.unit);  % 0.1 mm

nas = cfg.fiducial.nas;
ini = cfg.fiducial.ini;
lpa = cfg.fiducial.lpa;
rpa = cfg.fiducial.rpa;

if any(dist(headshape.pos, nas)<tolerance)
  ft_warning('Nasion coincides with headshape vertex, addding random displacement of about %f %s', tolerance, headshape.unit);
  nas = nas + tolerance*randn(1,3);
end
if any(dist(headshape.pos, ini)<tolerance)
  ft_warning('Inion coincides with headshape vertex, addding random displacement of about %f %s', tolerance, headshape.unit);
  ini = ini + tolerance*randn(1,3);
end
if any(dist(headshape.pos, lpa)<tolerance)
  ft_warning('LPA coincides with headshape vertex, addding random displacement of about %f %s', tolerance, headshape.unit);
  lpa = lpa + tolerance*randn(1,3);
end
if any(dist(headshape.pos, rpa)<tolerance)
  ft_warning('RPA coincides with headshape vertex, addding random displacement of about %f %s', tolerance, headshape.unit);
  rpa = rpa + tolerance*randn(1,3);
end

pos = headshape.pos;
tri = headshape.tri;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is largely shared with ELEC1020_LOCATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do not use the feedback of the low-level function
feedback = false;

% determine the approximate location of the vertex
ori = (lpa+rpa+nas+ini)/4;      % center of head
ver =  cross(rpa-lpa, nas-ini); % orientation
ver = ver /sqrt(norm(ver));     % make correct length
ver = ori + 0.7*ver;            % location from center of head

% point near LPA that is at 50% of left lower contour
[cnt1, cnt2] = elec1020_follow(pos, tri, nas, lpa, ini, feedback);
mle = elec1020_fraction(cnt1, cnt2, 0.5);

% point near RPA that is at 50% of right lower contour
[cnt1, cnt2] = elec1020_follow(pos, tri, nas, rpa, ini, feedback);
mre = elec1020_fraction(cnt1, cnt2, 0.5);

% determine two points that approximate the vertex
[cnt1, cnt2] = elec1020_follow(pos, tri, nas, ver, ini, feedback);
ver1 = elec1020_fraction(cnt1, cnt2, 0.5);
[cnt1, cnt2] = elec1020_follow(pos, tri, mle, ver, mre, feedback);
ver2 = elec1020_fraction(cnt1, cnt2, 0.5);

% refined estimate is the average of these two
ver = (ver1+ver2)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start contouring, note that the contour is constructed slightly higher
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ant-post contour through vertex
fprintf('constructing vertical ant-post contour\n');
[cnt1, cnt2] = elec1020_follow(pos, tri, nas, ver, ini, feedback);
Fpz  = elec1020_fraction(cnt1, cnt2,  2.5/20); %  2.5 instead of 2
Oz   = elec1020_fraction(cnt1, cnt2, 17.5/20); % 17.5 instead of 18

% left-right through vertex
fprintf('constructing C contour\n');
[cnt1, cnt2] = elec1020_follow(pos, tri, mle, ver, mre, feedback);
T7   = elec1020_fraction(cnt1, cnt2,  2.5/20);
T8   = elec1020_fraction(cnt1, cnt2, 17.5/20);

% horizontal ant-post through T7
fprintf('constructing horizontal left contour\n');
[left1, left2] = elec1020_follow(pos, tri, Fpz, T7, Oz, feedback);

% horizontal ant-post through T8
fprintf('constructing horizontal right contour\n');
[right1, right2] = elec1020_follow(pos, tri, Fpz, T8, Oz, feedback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project the left and right contour into a flat horozontal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nleft  = size(left1, 1);
nright = size(right1, 1);

all = [left1; left2; right1; right2];
meanpos = mean(all,1);
all(:,1) = all(:,1) - meanpos(1);
all(:,2) = all(:,2) - meanpos(2);
all(:,3) = all(:,3) - meanpos(3);

% flatten the contour
[u, s, v] = svd(all);
s(3,3) = 0;
all = u * s * v';
all(:,1) = all(:,1) + meanpos(1);
all(:,2) = all(:,2) + meanpos(2);
all(:,3) = all(:,3) + meanpos(3);

% split all points again in the four sets
left1  = all(1:(nleft),:);
left2  = all((nleft+1):(nleft+nleft),:);
right1 = all((nleft+nleft+1):(nleft+nleft+nright),:);
right2 = all((nleft+nleft+nright+1):(nleft+nleft+nright+nright),:);

% compute the length along the contours
circumference = 0;
for i=1:size(left1,1)
  circumference = circumference + dist(left1(i,:), left2(i,:));
end
for i=1:size(right1,1)
  circumference = circumference + dist(right1(i,:), right2(i,:));
end

% give graphical feedback
if istrue(cfg.feedback)
  figure
  ft_plot_mesh(struct('pos', pos, 'tri', tri), 'edgecolor', 'none', 'facecolor', 'skin')
  lighting gouraud
  material dull
  lightangle(0, 90);
  alpha 0.9
  grid on
  hold on
  view([1 1 0.5])

  % draw the anatomical landmarks
  ft_plot_mesh(nas, 'vertexsize', 30)
  ft_plot_mesh(lpa, 'vertexsize', 30)
  ft_plot_mesh(ini, 'vertexsize', 30)
  ft_plot_mesh(rpa, 'vertexsize', 30)
  ft_plot_mesh(ver, 'vertexsize', 30)

  % draw the contour over the left hemisphere
  X = [left1(:,1) left2(:,1)];
  Y = [left1(:,2) left2(:,2)];
  Z = [left1(:,3) left2(:,3)];
  line(X, Y, Z, 'Color', 'k', 'LineWidth', 3)

  % draw the contour over the right hemisphere
  X = [right1(:,1) right2(:,1)];
  Y = [right1(:,2) right2(:,2)];
  Z = [right1(:,3) right2(:,3)];
  line(X, Y, Z, 'Color', 'k', 'LineWidth', 3)
end

% give textual feedback
ft_info('the estimated head circumference is %f %s', circumference, headshape.unit);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
