function ft_plot_topo3d(pos, val, varargin)

% FT_PLOT_TOPO3D makes a 3-D topographic representation of the electric
% potential or field at the sensor locations
%
% Use as
%   ft_plot_topo3d(pos, val, ...);
% where the channel positions are given as a Nx3 matrix and the values are
% given as Nx1 vector.
%
% Optional input arguments should be specified in key-value pairs and can include
%   'contourstyle' = string, 'none', 'black', 'color' (default = 'none')
%   'isolines'     = vector with values at which to draw isocontours, or 'auto' (default = 'auto')
%   'facealpha'    = scalar, between 0 and 1 (default = 1)
%   'refine'       = scalar, number of refinement steps for the triangulation, to get a smoother interpolation (default = 0)
%
% See also FT_PLOT_TOPO, FT_TOPOPLOTER, FT_TOPOPLOTTFR

% Copyright (C) 2009-2015, Robert Oostenveld
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

ws = warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
contourstyle  = ft_getopt(varargin, 'contourstyle', 'none');
nrefine       = ft_getopt(varargin, 'refine', 0);
isolines      = ft_getopt(varargin, 'isolines', 'auto');
topostyle     = ft_getopt(varargin, 'topostyle', 'color');  % FIXME what is the purpose of this option?
facealpha     = ft_getopt(varargin, 'facealpha', 1);

if islogical(contourstyle) && contourstyle==false
  % false was supported up to 18 November 2013, 'none' is more consistent with other plotting options
  contourstyle = 'none';
end

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

if size(val,2)==size(pos,1)
  val = val';
end

% the interpolation requires a triangulation
tri = projecttri(pos, 'delaunay');

if nrefine>0,
  posorig = pos;
  triorig = tri;
  valorig = val;
  for k = 1:nrefine
    [pos,tri] = refine(pos, tri);
  end
  prjorig = elproj(posorig);
  prj     = elproj(pos);
  val     = griddata(prjorig(:,1),prjorig(:,2),valorig,prj(:,1),prj(:,2),'v4');
  if numel(facealpha)==size(posorig,1)
    facealpha = griddata(prjorig(:,1),prjorig(:,2),facealpha,prj(:,1),prj(:,2),'v4');
  end
end


if ~isequal(topostyle, false)
  switch topostyle
    case 'color'
      % plot a 2D or 3D triangulated surface with linear interpolation
      if length(val)==size(pos,1)
        hs = patch('Vertices', pos, 'Faces', tri, 'FaceVertexCData', val, 'FaceColor', 'interp');
      else
        hs = patch('Vertices', pos, 'Faces', tri, 'CData', val, 'FaceColor', 'flat');
      end
      set(hs, 'EdgeColor', 'none');
      set(hs, 'FaceLighting', 'none');
      
      % if facealpha is an array with number of elements equal to the number of vertices
      if size(pos,1)==numel(facealpha)
        set(hs, 'FaceVertexAlphaData', facealpha);
        set(hs, 'FaceAlpha', 'interp');
      elseif ~isempty(pos) && numel(facealpha)==1 && facealpha~=1
        % the default is 1, so that does not have to be set
        set(hs, 'FaceAlpha', facealpha);
      end
      
    otherwise
      error('unsupported topostyle');
  end % switch contourstyle
end % plot the interpolated topography


if ~strcmp(contourstyle, 'none')
  
  if ischar(isolines)
    if isequal(isolines, 'auto')
      minval = min(val);
      maxval = max(val);
      scale = max(abs(minval), abs(maxval));
      scale = 10^(floor(log10(scale))-1);
      minval = floor(minval/scale)*scale;
      maxval = ceil(maxval/scale)*scale;
      isolines = minval:scale:maxval;
    else
      error('unsupported isolines');
    end
  end % convert string to vector
  
  tri_val = val(tri);
  tri_min = min(tri_val, [], 2);
  tri_max = max(tri_val, [], 2);
  
  for cnt_indx=1:length(isolines)
    cnt = isolines(cnt_indx);
    use = cnt>=tri_min & cnt<=tri_max;
    counter = 0;
    intersect1 = [];
    intersect2 = [];
    
    for tri_indx=find(use)'
      tri_pos = pos(tri(tri_indx,:), :);
      v(1) = tri_val(tri_indx,1);
      v(2) = tri_val(tri_indx,2);
      v(3) = tri_val(tri_indx,3);
      la(1) = (cnt-v(1)) / (v(2)-v(1)); % abcissa between vertex 1 and 2
      la(2) = (cnt-v(2)) / (v(3)-v(2)); % abcissa between vertex 2 and 3
      la(3) = (cnt-v(3)) / (v(1)-v(3)); % abcissa between vertex 1 and 2
      abc(1,:) = tri_pos(1,:) + la(1) * (tri_pos(2,:) - tri_pos(1,:));
      abc(2,:) = tri_pos(2,:) + la(2) * (tri_pos(3,:) - tri_pos(2,:));
      abc(3,:) = tri_pos(3,:) + la(3) * (tri_pos(1,:) - tri_pos(3,:));
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
  
  % collect all different contour isolines for plotting
  intersect1 = [];
  intersect2 = [];
  cntlevel   = [];
  for cnt_indx=1:length(isolines)
    intersect1 = [intersect1; contour(cnt_indx).intersect1];
    intersect2 = [intersect2; contour(cnt_indx).intersect2];
    cntlevel   = [cntlevel; ones(contour(cnt_indx).n,1) * isolines(cnt_indx)];
  end
  
  X = [intersect1(:,1) intersect2(:,1)]';
  Y = [intersect1(:,2) intersect2(:,2)]';
  C = [cntlevel(:)     cntlevel(:)]';
  
  if size(pos,2)>2
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

warning(ws); % revert to original state
