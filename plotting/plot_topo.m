function [varargout] = plot_topo(chanX, chanY, dat, varargin)

% PLOT_TOPO interpolates and plots the 2-D spatial topography of the
% potential or field distribution over the head
%
% Use as
%   plot_topo(x, y, val, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'hpos'
%   'vpos'
%   'width'
%   'height'
%   'shading'
%   'gridscale'
%   'mask'
%   'outline'

% Copyrights (C) 2009, Giovanni Piantoni
%
% $Log: plot_topo.m,v $
% Revision 1.9  2009/10/09 11:18:23  jansch
% fixed some inappropriate behaviour
%
% Revision 1.8  2009/10/09 10:23:34  jansch
% added interplim as option, according to topoplot (default = electrodes)
%
% Revision 1.7  2009/08/12 15:15:18  jansch
% also changed the order of inputs on the first line of the function
%
% Revision 1.6  2009/08/05 08:58:54  roboos
% changed the order of the input arguments to plot_topo from (val, x, y) into (x, y, val)
%
% Revision 1.5  2009/08/05 08:53:20  roboos
% plot the outline of the head if specified
% keep hold on/off the same
%
% Revision 1.4  2009/07/29 15:04:16  giopia
% resolved ambiguity of var mask
%
% Revision 1.3  2009/07/29 10:24:24  roboos
% construct the binary image for masking inside this function and reuse it as long as the relevant input does not change
% this is achieved with a persistent variable and by checking the input arguments
%
% Revision 1.2  2009/06/02 15:36:25  giopia
% first implementation based on topoplot.m
%

% these are for speeding up the plotting on subsequent calls
persistent previous_argin previous_maskimage

holdflag = ishold;
hold on

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'gridscale', 'shading', 'mask', 'outline', 'interplim'});
hpos        = keyval('hpos',      varargin); if isempty(hpos);       hpos = 0;           end
vpos        = keyval('vpos',      varargin); if isempty(vpos);       vpos = 0;           end
width       = keyval('width',     varargin); if isempty(width);      width = 1;          end
height      = keyval('height',    varargin); if isempty(height);     height = 1;         end
gridscale   = keyval('gridscale', varargin); if isempty(gridscale);  gridscale = 67;     end; % 67 in original
shading     = keyval('shading',   varargin); if isempty(shading);    shading = 'flat';   end;
mask        = keyval('mask',      varargin);
outline     = keyval('outline',   varargin);
interplim   = keyval('interplim', varargin); if isempty(interplim);  interplim = 'electrodes'; end

chanX = chanX * width  + hpos;
chanY = chanY * height + vpos;

if strcmp(interplim, 'electrodes'),
  hlim = [min(chanX) max(chanX)];
  vlim = [min(chanY) max(chanY)];
elseif strcmp(interplim, 'mask') && ~isempty(mask),
  hlim = [inf -inf];
  vlim = [inf -inf];
  for i=1:length(mask)
    hlim = [min([hlim(1); mask{i}(:,1)+hpos]) max([hlim(2); mask{i}(:,1)+hpos])];
    vlim = [min([vlim(1); mask{i}(:,2)+vpos]) max([vlim(2); mask{i}(:,2)+vpos])];
  end
else
  hlim = [min(chanX) max(chanX)];
  vlim = [min(chanY) max(chanY)];
end

% try to speed up the preparation of the mask on subsequent calls
current_argin = {chanX, chanY, gridscale, mask};
if isequal(current_argin, previous_argin)
  % don't construct the binary image, but reuse it from the previous call
  maskimage = previous_maskimage;
elseif ~isempty(mask)
  % convert the mask into a binary image
  maskimage = false(gridscale);
  %hlim      = [min(chanX) max(chanX)];
  %vlim      = [min(chanY) max(chanY)];
  xi        = linspace(hlim(1), hlim(2), gridscale);   % x-axis for interpolation (row vector)
  yi        = linspace(vlim(1), vlim(2), gridscale);   % y-axis for interpolation (row vector)
  [Xi,Yi]   = meshgrid(xi', yi);
  for i=1:length(mask)
    mask{i}(:,1) = mask{i}(:,1)+hpos;
    mask{i}(:,2) = mask{i}(:,2)+vpos;
    mask{i}(end+1,:) = mask{i}(1,:);                   % force them to be closed
    maskimage(inside_contour([Xi(:) Yi(:)], mask{i})) = true;
  end
else
  maskimage = [];
end

xi         = linspace(hlim(1), hlim(2), gridscale);       % x-axis for interpolation (row vector)
yi         = linspace(vlim(1), vlim(2), gridscale);       % y-axis for interpolation (row vector)
[Xi,Yi,Zi] = griddata(chanX', chanY, dat, xi', yi, 'v4'); % interpolate the topographic data

if ~isempty(maskimage)
  % apply anatomical mask to the data, i.e. that determines that the interpolated data outside the circle is not displayed
  Zi(~maskimage) = NaN;
end

deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry
h = surface(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)), Zi, 'EdgeColor', 'none', 'FaceColor', shading);

% plot the outline of the head, ears and nose
for i=1:length(outline)
  xval = outline{i}(:,1) * width  + hpos;
  yval = outline{i}(:,2) * height + vpos;
  plot(xval, yval, 'Color', 'k', 'LineWidth', 1)
end

% the (optional) output is the handle
if nargout == 1
  varargout{1} = h;
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
if isempty(previous_argin)
  previous_argin     = current_argin;
  previous_maskimage = maskimage;
end

if ~holdflag
  hold off
end

