function ft_plot_topo3d(pnt, val, varargin)

% FT_PLOT_TOPO3D makes a 3-D topographic representation of the electric
% potential or field at the sensor locations
%
% Use as
%   ft_plot_topo3d(pos, val, ...);
% where the channel positions are given as a Nx3 matrix and the values are
% given as Nx1 vector.
%
% Optional input arguments should be specified in key-value pairs and can include
% 'contourstyle'  false  (default), 'black', 'color' makes contours of b/w or colored lines
% 'isocontour'    'auto' (default - and only - option)
% 'topostyle'     'color'(default - and only - option)

%
% See also FT_PLOT_TOPO
%
% Copyright (C) 2009, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

ws = warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
contourstyle  = ft_getopt(varargin, 'contourstyle', false);
isocontour    = ft_getopt(varargin, 'isocontour',   'auto');
topostyle     = ft_getopt(varargin, 'topostyle',    'color');

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

% the interpolation requires a triangulation
tri = projecttri(pnt, 'delaunay');

if ~isequal(topostyle, false)
  switch topostyle
    case 'color'
      % plot a 2D or 3D triangulated surface with linear interpolation
      if length(val)==size(pnt,1)
        hs = patch('Vertices', pnt, 'Faces', tri, 'FaceVertexCData', val, 'FaceColor', 'interp');
      else
        hs = patch('Vertices', pnt, 'Faces', tri, 'CData', val, 'FaceColor', 'flat');
      end
      set(hs, 'EdgeColor', 'none');
      set(hs, 'FaceLighting', 'none');
    otherwise
      error('unsupported topostyle');
  end % switch contourstyle
end % plot the interpolated topography


if ~isequal(contourstyle, false)
  
  if isequal(isocontour, 'auto')
    minval = min(val);
    maxval = max(val);
    scale = max(abs(minval), abs(maxval));
    scale = 10^(floor(log10(scale))-1);
    minval = floor(minval/scale)*scale;
    maxval = ceil(maxval/scale)*scale;
    isocontour = minval:scale:maxval;
  end
  
  triangle_val = val(tri);
  triangle_min = min(triangle_val, [], 2);
  triangle_max = max(triangle_val, [], 2);
  
  for cnt_indx=1:length(isocontour)
    cnt = isocontour(cnt_indx);
    use = cnt>=triangle_min & cnt<=triangle_max;
    counter = 0;
    intersect1 = [];
    intersect2 = [];
    
    for tri_indx=find(use)'
      pos  = pnt(tri(tri_indx,:), :);
      v(1) = triangle_val(tri_indx,1);
      v(2) = triangle_val(tri_indx,2);
      v(3) = triangle_val(tri_indx,3);
      la(1) = (cnt-v(1)) / (v(2)-v(1)); % abcissa between vertex 1 and 2
      la(2) = (cnt-v(2)) / (v(3)-v(2)); % abcissa between vertex 2 and 3
      la(3) = (cnt-v(3)) / (v(1)-v(3)); % abcissa between vertex 1 and 2
      abc(1,:) = pos(1,:) + la(1) * (pos(2,:) - pos(1,:));
      abc(2,:) = pos(2,:) + la(2) * (pos(3,:) - pos(2,:));
      abc(3,:) = pos(3,:) + la(3) * (pos(1,:) - pos(3,:));
      counter = counter + 1;
      sel     = find(la>=0 & la<=1);
      intersect1(counter, :) = abc(sel(1),:);
      intersect2(counter, :) = abc(sel(2),:);
    end
    
    % remember the details for external reference
    contour(cnt_indx).level = cnt;
    contour(cnt_indx).n     = counter;
    contour(cnt_indx).intersect1 = intersect1;
    contour(cnt_indx).intersect2 = intersect2;
  end
  
  % collect all different contour isocontour for plotting
  intersect1 = [];
  intersect2 = [];
  cntlevel   = [];
  for cnt_indx=1:length(isocontour)
    intersect1 = [intersect1; contour(cnt_indx).intersect1];
    intersect2 = [intersect2; contour(cnt_indx).intersect2];
    cntlevel   = [cntlevel; ones(contour(cnt_indx).n,1) * isocontour(cnt_indx)];
  end
  
  X = [intersect1(:,1) intersect2(:,1)]';
  Y = [intersect1(:,2) intersect2(:,2)]';
  C = [cntlevel(:)     cntlevel(:)]';
  
  if size(pnt,2)>2
    Z = [intersect1(:,3) intersect2(:,3)]';
  else
    Z = zeros(2, length(cntlevel));
  end
  
  switch contourstyle
    case 'black'
      % make black-white contours
      hc = [];
      for i=1:length(cntlevel)
        if cntlevel(i)>0
          linestyle = '-';
          linewidth = 1;
        elseif cntlevel(i)<0
          linestyle = '--';
          linewidth = 1;
        else
          linestyle = '-';
          linewidth = 2;
        end
        h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
          'ZData', Z(:,i), 'CData', C(:,i), ...
          'facecolor','none','edgecolor','black', ...
          'linestyle', linestyle, 'linewidth', linewidth, ...
          'userdata',cntlevel(i));
        hc = [hc; h1];
      end
      
    case 'color'
      % make full-color contours
      hc = [];
      for i=1:length(cntlevel)
        h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
          'ZData', Z(:,i), 'CData', C(:,i), ...
          'facecolor','none','edgecolor','flat',...
          'userdata',cntlevel(i));
        hc = [hc; h1];
      end
      
    otherwise
      error('unsupported contourstyle');
  end % switch contourstyle
  
end % plot the contours

axis off
axis vis3d
axis equal

if ~holdflag
  hold off
end

warning(ws); %revert to original state
