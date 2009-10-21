function [h, lpa, rpa] = volplot(x, y, z, dat, sel, cscale)

% VOLPLOT make 2D or 3D plot of volumetric data (e.g. MRI)
% that is defined on a regular orthogonal grid
% 
% volplot(dat, sel) or
% volplot(x, y, z, dat, sel)
% volplot(x, y, z, dat, sel, caxis)
%
% where sel is one of
%   [x, y, z]     intersection through the three orthogonal directions
%   index         linear index of the voxel of interest
%   'min'         intersection at the minimum
%   'max'         intersection at the maximum
%   'center'      intersect at the center of each axis
%   'interactive' intersect at the center, then go into interactive mode
%   'maxproject'  project the maximum value along each orthogonal direction
%   'sumproject'  integrated value along each orthogonal direction (glassbrain)
%   'montage'     show all slices
% and caxis is the [min max] used for the color scaling
% 
% See also TRIPLOT, LINEPLOT (in ~roberto/matlab/misc)
% See also NDGRID

% Copyright (C) 2003, Robert Oostenveld
% 
% $Log: volplot.m,v $
% Revision 1.15  2007/01/04 12:26:11  roboos
% added linear index of the voxel as selection method
%
% Revision 1.14  2006/05/02 19:14:19  roboos
% default to interactive when only one 3D input
%
% Revision 1.13  2004/06/28 07:52:42  roberto
% fixed some minor bugs, changed crosshair color from black to yello
%
% Revision 1.12  2004/01/29 09:14:03  roberto
% added an interactive mode to localize nas/lpa/rpa
%
% Revision 1.11  2004/01/21 12:55:49  roberto
% added voxel number to printed output, easy find back maximum source
%
% Revision 1.10  2004/01/19 14:50:13  roberto
% added option 'center'
%
% Revision 1.9  2004/01/13 09:24:13  roberto
% added input option for color axis
%
% Revision 1.8  2003/11/03 11:42:36  roberto
% added SUMPROJECT (i.e. glassbrain type of projection)
% fixed bug with fprintf of uint8
%
% Revision 1.7  2003/09/22 10:33:39  roberto
% added fprintf statement indicating value and position
%
% Revision 1.6  2003/09/03 09:53:21  roberto
% removed automatic scaling of montage display
%
% Revision 1.5  2003/04/23 10:20:49  roberto
% fixed 2 bugs related to 1D input data
%
% Revision 1.4  2003/03/24 13:15:04  roberto
% added support for 1D data (automatically converted to 3D)
%
% Revision 1.3  2003/03/21 14:04:31  roberto
% added axis xy to montage plot
%
% Revision 1.2  2003/03/17 10:23:04  roberto
% added montage display
% added automatic min/max selection of intersection point
%

if nargin<2
  dat = x;
  sel = 'interactive'; 
  x = 1:size(dat,1);
  y = 1:size(dat,2);
  z = 1:size(dat,3);
elseif nargin<3
  dat = x;
  sel = y;
  x = 1:size(dat,1);
  y = 1:size(dat,2);
  z = 1:size(dat,3);
elseif nargin>4
  if length(size(dat))==3
    if length(x)~=size(dat,1), error('incorrect x-axis specification'), end
    if length(y)~=size(dat,2), error('incorrect y-axis specification'), end
    if length(z)~=size(dat,3), error('incorrect z-axis specification'), end
  end
else
  error('incorrect number of input arguments');
end

if nargin==6
  cmin = cscale(1);
  cmax = cscale(2);
else
  % ensure same color scaling for all figures
  cmin = min(dat(:));
  cmax = max(dat(:));
end
  
% convert from 1-d to 3-d array
if any(size(dat)==1)
  dim(1) = length(x);
  dim(2) = length(y);
  dim(3) = length(z);
  dat = reshape(dat, dim);
else
  dim = size(dat);
end

% convert the selection to the indices of the x/y/z intersection 
if isstr(sel) & strcmp(sel, 'min')
  [minval, minindx] = min(dat(:));
  [xi, yi, zi] = ind2sub(size(dat), minindx);
elseif isstr(sel) & strcmp(sel, 'max')
  [maxval, maxindx] = max(dat(:));
  [xi, yi, zi] = ind2sub(size(dat), maxindx);
elseif isstr(sel) & strcmp(sel, 'center')
  xi = round(length(x)/2);
  yi = round(length(y)/2);
  zi = round(length(z)/2);
elseif isstr(sel) & strcmp(sel, 'interactive')
  xi = round(length(x)/2);
  yi = round(length(y)/2);
  zi = round(length(z)/2);
elseif ~isstr(sel) && length(sel)==1
  [xi, yi, zi] = ind2sub(dim, sel);
else
  xi = nearest(x, sel(1));
  yi = nearest(y, sel(2));
  zi = nearest(z, sel(3));
end

% start the plotting
if strcmp(sel, 'interactive')
  xc = x(xi);
  yc = y(yi);
  zc = z(zi);
  nas = [];
  lpa = [];
  rpa = [];
  while(1)
    fprintf('============================================================\n');
    fprintf('click with mouse button to reslice the display to a new position\n');
    fprintf('press n/l/r on keyboard to record a fiducial position\n');
    fprintf('press q on keyboard to quit interactive mode\n');
    volplot(x, y, z, dat, [xc yc zc]);
    drawnow;
    try, [d1, d2, key] = ginput(1); catch, key='q'; end
    if key=='q'
      % rename the first output argument of this function
      h = nas;
      break;
    elseif key=='l'
      lpa = [xc yc zc];
    elseif key=='r'
      rpa = [xc yc zc];
    elseif key=='n'
      nas = [xc yc zc];
    else         
      % update the view to a new position
      l1 = get(get(gca, 'xlabel'), 'string');
      l2 = get(get(gca, 'ylabel'), 'string');
      switch l1,
        case 'x'
          xc = d1;
        case 'y'
          yc = d1;
        case 'z'
          zc = d1;
      end
      switch l2,
        case 'x'
          xc = d2;
        case 'y'
          yc = d2;
        case 'z'
          zc = d2;
      end
    end
    if ~isempty(nas), fprintf('nas = [%f %f %f]\n', nas); else fprintf('nas = undefined\n'); end 
    if ~isempty(lpa), fprintf('lpa = [%f %f %f]\n', lpa); else fprintf('lpa = undefined\n'); end 
    if ~isempty(rpa), fprintf('rpa = [%f %f %f]\n', rpa); else fprintf('rpa = undefined\n'); end 
  end
  
elseif strcmp(sel, 'montage')
  % make plot of x-y slices for all z values
  maxval = max(dat(:));
  for z=1:size(dat,3)
    % convert to 4D image for montage display
    % transpose to correct for x-y axis change in Matlab image function
    img(:,:,1,z) = transpose(dat(:,:,z));
  end
  montage(img);
  axis xy

elseif strcmp(sel, 'sumproject')
  % make plot of integrated-value projection along the thee orthogonal directions
  delete(subplot(2,2,4));	% delete the old colorbar
  h1 = subplot(2,2,1);
  h2 = subplot(2,2,2);
  h3 = subplot(2,2,3);
  h4 = subplot(2,2,4);	% this will  be the new colorbar

  % change not-a-number values to zero
  dat(find(isnan(dat(:)))) = 0;

  subplot(h1);
  imagesc(x, z, squeeze(sum(dat, 2))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('z');
  caxis([cmin cmax]);

  subplot(h2);
  imagesc(y, z, squeeze(sum(dat, 1))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('y'); ylabel('z');
  caxis([cmin cmax]);

  subplot(h3);
  imagesc(x, y, squeeze(sum(dat, 3))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('y');
  caxis([cmin cmax]);

  subplot(h4);
  colorbar(h4, 'peer', h1);
  xlabel('colorscale')

elseif strcmp(sel, 'maxproject')
  % make plot of maximum-value projection along the thee orthogonal directions
  delete(subplot(2,2,4));	% delete the old colorbar
  h1 = subplot(2,2,1);
  h2 = subplot(2,2,2);
  h3 = subplot(2,2,3);
  h4 = subplot(2,2,4);	% this will  be the new colorbar

  subplot(h1);
  imagesc(x, z, squeeze(max(dat, [], 2))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('z');
  caxis([cmin cmax]);

  subplot(h2);
  imagesc(y, z, squeeze(max(dat, [], 1))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('y'); ylabel('z');
  caxis([cmin cmax]);

  subplot(h3);
  imagesc(x, y, squeeze(max(dat, [], 3))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('y');
  caxis([cmin cmax]);

  subplot(h4);
  colorbar(h4, 'peer', h1);
  xlabel('colorscale')

else
  % make plot of three orthogonal slices intersecting at [xi yi zi]
  if ~exist('xi', 'var') | ~exist('yi', 'var') | ~exist('zi', 'var')
    error('nothing to plot, no selection given')
  end

  fprintf('value of %f in voxel %d at [%.02f %.02f %.02f]\n', double(dat(xi, yi, zi)), sub2ind(dim, xi, yi, zi), x(xi), y(yi), z(zi));

  delete(subplot(2,2,4));	% delete the old colorbar
  h1 = subplot(2,2,1);
  h2 = subplot(2,2,2);
  h3 = subplot(2,2,3);
  h4 = subplot(2,2,4);	% this will  be the new colorbar

  subplot(h1);
  imagesc(x, z, squeeze(dat(:,yi,:))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('z');
  caxis([cmin cmax]);
  crosshair([x(xi) z(zi)], 'color', 'yellow');

  subplot(h2);
  imagesc(y, z, squeeze(dat(xi,:,:))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('y'); ylabel('z');
  caxis([cmin cmax]);
  crosshair([y(yi) z(zi)], 'color', 'yellow');

  subplot(h3);
  imagesc(x, y, squeeze(dat(:,:,zi))'); set(gca, 'ydir', 'normal')
  axis equal; axis tight;
  xlabel('x'); ylabel('y');
  caxis([cmin cmax]);
  crosshair([x(xi) y(yi)], 'color', 'yellow');

  subplot(h4);
  colorbar(h4, 'peer', h1);
  xlabel('colorscale')
end

