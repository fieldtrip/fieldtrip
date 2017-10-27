function ft_plot_montage(dat, varargin)

% FT_PLOT_MONTAGE makes a montage of a 3-D array by selecting slices at regular distances
% and combining them in one large 2-D image.  Note that the montage of MRI slices is not to
% be confused with the EEG montage, which is a way of specifying the reference scheme
% between electrodes.
%
% Use as
%   ft_plot_montage(dat, ...)
% where dat is a 3-D array.
% 
% Additional options should be specified in key-value pairs and can be
%   'transform'     = 4x4 homogeneous transformation matrix specifying the mapping from voxel space to the coordinate system in which the data are plotted.
%   'location'      = 1x3 vector specifying a point on the plane which will be plotted the coordinates are expressed in the coordinate system in which the data will be plotted. location defines the origin of the plane
%   'orientation'   = 1x3 vector specifying the direction orthogonal through the plane which will be plotted (default = [0 0 1])
%   'srange'        = 
%   'slicesize'     = 
%   'nslice'        = scalar, number of slices
%   'maskstyle'     = string, 'opacity' or 'colormix', defines the rendering
%   'background'    = needed when maskstyle is 'colormix', 3D-matrix with
%                     the same size as the data matrix, serving as
%                     grayscale image that provides the background
% 
% See also FT_PLOT_ORTHO, FT_PLOT_SLICE, FT_SOURCEPLOT

% undocumented, these are passed on to FT_PLOT_SLICE
%   'intersectmesh'  = triangulated mesh through which the intersection of the plane will be plotted (e.g. cortical sheet)
%   'intersectcolor' = color for the intersection

% Copyrights (C) 2012, Jan-Mathijs Schoffelen
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

transform       = ft_getopt(varargin, 'transform', eye(4));
loc             = ft_getopt(varargin, 'location');
ori             = ft_getopt(varargin, 'orientation');
srange          = ft_getopt(varargin, 'slicerange');
slicesize       = ft_getopt(varargin, 'slicesize');
nslice          = ft_getopt(varargin, 'nslice');
backgroundcolor = ft_getopt(varargin, 'backgroundcolor', [0 0 0]);

% the intersectmesh and plotmarker options are passed on to FT_PLOT_SLICE
dointersect = ~isempty(ft_getopt(varargin, 'intersectmesh'));
domarker    = ~isempty(ft_getopt(varargin, 'plotmarker'));

% set the location if empty
if isempty(loc) && (isempty(transform) || isequal(transform, eye(4)))
  % go to the middle of the volume if the data seem to be in voxel coordinates
  loc = size(dat)./2;
elseif isempty(loc)
  % otherwise take the origin of the coordinate system
  loc = [0 0 0];
end

% check compatibility of inputs
if size(loc, 1) == 1 && isempty(nslice)
  nslice = 20;
elseif size(loc, 1) == 1 && ~isempty(nslice)
  % this is not a problem, slice spacing will be determined
elseif size(loc, 1) > 1 && isempty(nslice)
  % this is not a problem, number of slices is determined by loc
  nslice = size(loc, 1);
elseif size(loc, 1) > 1 && ~isempty(nslice)
  if size(loc, 1) ~= nslice
    ft_error('you should either specify a set of locations or a single location with a number of slices');
  end
end

% set the orientation if empty
if isempty(ori)
  ori = [0 0 1];
end

% ensure the ori to have unit norm
for k = 1:size(ori,1)
  ori(k,:) = ori(k,:)./norm(ori(k,:));
end

% determine the slice range
if size(loc, 1) == 1 && nslice > 1
  if isempty(srange) || (ischar(srange) && strcmp(srange, 'auto'))
    srange = [-50 70];
  else
  end
  loc = repmat(loc, [nslice 1]) + linspace(srange(1),srange(2),nslice)'*ori;
end

% ensure that the ori has the same size as the loc
if size(ori,1)==1 && size(loc,1)>1
  ori = repmat(ori, size(loc,1), 1);
end

div     = [ceil(sqrt(nslice)) ceil(sqrt(nslice))];
optarg  = varargin;
corners = [inf -inf inf -inf inf -inf]; % get the corners for the axis specification

for k = 1:nslice
  % define 'x' and 'y' axis in projection plane, the definition of x and y is more or less arbitrary
  [x, y] = projplane(ori(k,:)); % z = ori
  
  % get the transformation matrix to project onto the xy-plane
  T  = [x(:) y(:) ori(k,:)' loc(k,:)'; 0 0 0 1];
  
  optarg = ft_setopt(optarg, 'location',    loc(k,:));
  optarg = ft_setopt(optarg, 'orientation', ori(k,:));
  ix     = mod(k-1, div(1));
  iy     = floor((k-1)/div(1));
  h(k)   = ft_plot_slice(dat, optarg{:}); % FIXME is it safe to pass all optinoal inputs?
  
  xtmp = get(h(k), 'xdata');
  ytmp = get(h(k), 'ydata');
  ztmp = get(h(k), 'zdata');
  siz  = size(xtmp);
  if k==1 && isempty(slicesize)
    slicesize = siz;
  end
  
  % project the positions onto the xy-plane
  pos = [xtmp(:) ytmp(:) ztmp(:)];
  pos = ft_warp_apply(inv(T), pos);
  
  xtmp = reshape(pos(:,1), siz);
  ytmp = reshape(pos(:,2), siz);
  ztmp = reshape(pos(:,3), siz);
  
  % add some offset in the x and y directions to create the montage
  offset(1) = iy*(slicesize(1)-1);
  offset(2) = ix*(slicesize(2)-1); 
  
  % update the specification of the corners of the montage plot
  if ~isempty(xtmp)
    c1 = offset(1) + min(xtmp(:));
    c2 = offset(1) + max(xtmp(:));
    c3 = offset(2) + min(ytmp(:));
    c4 = offset(2) + max(ytmp(:));
    c5 = min(ztmp(:));
    c6 = max(ztmp(:));
  end
  corners = [min(corners(1),c1) max(corners(2),c2) min(corners(3),c3) max(corners(4),c4) min(corners(5),c5) max(corners(6),c6)];
  
  % update the positions
  set(h(k), 'ydata', offset(1) + xtmp);
  set(h(k), 'xdata', offset(2) + ytmp);
  set(h(k), 'zdata',         0 * ztmp);
  
  if dointersect || domarker
    if ~exist('pprevious', 'var'), pprevious = []; end
    p = setdiff(findobj(gcf, 'type', 'patch'), pprevious);
    for kk = 1:numel(p)
      xtmp = get(p(kk), 'xdata');
      ytmp = get(p(kk), 'ydata');
      ztmp = get(p(kk), 'zdata');
      siz2 = size(xtmp);
      
      pos = [xtmp(:) ytmp(:) ztmp(:)];
      pos = ft_warp_apply(inv(T), pos);
  
      xtmp = reshape(pos(:,1), siz2);
      ytmp = reshape(pos(:,2), siz2);
      ztmp = reshape(pos(:,3), siz2);
  
      % update the positions
      set(p(kk), 'ydata', offset(1) + xtmp);
      set(p(kk), 'xdata', offset(2) + ytmp);
      set(p(kk), 'zdata',         0.0001 * ztmp);
    end
    pprevious = [pprevious(:);p(:)];
  end
  % drawnow; %this statement slows down the process big time on some file
  %systems. I don't know what's going on there, but the statement is not
  %really necessary, so commented out.
end
set(gcf, 'color', backgroundcolor);
set(gca, 'zlim', [0 1]);
%axis equal;
axis off;
view([0 90]);
axis(corners([3 4 1 2]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = projplane(z)
[u, s, v] = svd([eye(3) z(:)]);
x = u(:, 2)';
y = u(:, 3)';
