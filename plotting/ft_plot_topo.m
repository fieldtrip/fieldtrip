function Zi = ft_plot_topo(chanX, chanY, dat, varargin)

% FT_PLOT_TOPO interpolates and plots the 2-D spatial topography of the
% potential or field distribution over the head
%
% Use as
%   ft_plot_topo(x, y, val, ...)
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
%   'isolines'
%   'interplim'
%   'interpmethod'
%   'style'
%   'datmask'

% Copyrights (C) 2009, Giovanni Piantoni
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

% these are for speeding up the plotting on subsequent calls
persistent previous_argin previous_maskimage

warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'gridscale', 'shading', 'mask', 'outline', 'interplim', 'interpmethod','isolines','style', 'datmask'});
hpos          = keyval('hpos',         varargin);    if isempty(hpos);         hpos = 0;                 end
vpos          = keyval('vpos',         varargin);    if isempty(vpos);         vpos = 0;                 end
width         = keyval('width',        varargin);    if isempty(width);        width = 1;                end
height        = keyval('height',       varargin);    if isempty(height);       height = 1;               end
gridscale     = keyval('gridscale',    varargin);    if isempty(gridscale);    gridscale = 67;           end; % 67 in original
shading       = keyval('shading',      varargin);    if isempty(shading);      shading = 'flat';         end;
mask          = keyval('mask',         varargin);
outline       = keyval('outline',      varargin);
interplim     = keyval('interplim',    varargin);    if isempty(interplim);    interplim = 'electrodes'; end
interpmethod  = keyval('interpmethod', varargin);    if isempty(interpmethod); interpmethod = 'v4';      end
isolines      = keyval('isolines',     varargin);      
style         = keyval('style',        varargin);    if isempty(style);        style = 'surfiso';       end % can be 'surf', 'iso', 'isofill', 'surfiso'
datmask       = keyval('datmask',      varargin);

% everything is added to the current figure
holdflag = ishold;
hold on

chanX = chanX(:) * width  + hpos;
chanY = chanY(:) * height + vpos;

if strcmp(interplim, 'electrodes'),
  hlim = [min(chanX) max(chanX)];
  vlim = [min(chanY) max(chanY)];
elseif strcmp(interplim, 'mask') && ~isempty(mask),
  hlim = [inf -inf];
  vlim = [inf -inf];
  for i=1:length(mask)
    hlim = [min([hlim(1); mask{i}(:,1)*width+hpos]) max([hlim(2); mask{i}(:,1)*width+hpos])];
    vlim = [min([vlim(1); mask{i}(:,2)*width+vpos]) max([vlim(2); mask{i}(:,2)*width+vpos])];
  end
else
  hlim = [min(chanX) max(chanX)];
  vlim = [min(chanY) max(chanY)];
end

% check if all mask point are inside the limits otherwise redefine mask
newpoints = [];
if length(mask)==1
  % which channels are outside
  outside = false(length(chanX),1);
  inside  = inside_contour([chanX chanY], mask{1});
  outside = ~inside;
  newpoints = [chanX(outside) chanY(outside)];
end


% try to speed up the preparation of the mask on subsequent calls
current_argin = {chanX, chanY, gridscale, mask, datmask};
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
  if ~isempty(newpoints) && (hpos == 0 || vpos == 0)
    warning('Some points fall outside the outline, please consider using another layout')
% FIXME: I am not sure about it, to be tested!
%     tmp = [mask{1};newpoints];
%     indx = convhull(tmp(:,1),tmp(:,2));
%     mask{1} = tmp(indx,:);
% NOTE: if you set hpos and/or vpos, newpoints is not empty, but nothing
% needs to be fixed (this fixme screws up things, then)
  end 
  for i=1:length(mask)
    mask{i}(:,1) = mask{i}(:,1)*width+hpos;
    mask{i}(:,2) = mask{i}(:,2)*height+vpos;
    mask{i}(end+1,:) = mask{i}(1,:);                   % force them to be closed
    maskimage(inside_contour([Xi(:) Yi(:)], mask{i})) = true;
  end

else
  maskimage = [];
end

% adjust maskimage to also mask channels as specified in maskdat
if ~isempty(datmask)
  xi           = linspace(hlim(1), hlim(2), gridscale);   % x-axis for interpolation (row vector)
  yi           = linspace(vlim(1), vlim(2), gridscale);   % y-axis for interpolation (row vector)
  maskimagetmp = griddata(chanX', chanY, datmask, xi', yi, interpmethod); % interpolate the mask data
  if isempty(maskimage)
    maskimage = maskimagetmp;
  else
    maskimagetmp2 = maskimage + maskimagetmp;
    maskimage = maskimagetmp2 > 1;
  end
end
  
xi         = linspace(hlim(1), hlim(2), gridscale);       % x-axis for interpolation (row vector)
yi         = linspace(vlim(1), vlim(2), gridscale);       % y-axis for interpolation (row vector)
[Xi,Yi,Zi] = griddata(chanX', chanY, dat, xi', yi, interpmethod); % interpolate the topographic data

if ~isempty(maskimage)
  % apply anatomical mask to the data, i.e. that determines that the interpolated data outside the circle is not displayed
  Zi(~maskimage) = NaN;
end

if exist('maskimagetmp')
  maskimagetmp(~maskimage) = NaN;
end

% plot the outline of the head, ears and nose
for i=1:length(outline)
  xval = outline{i}(:,1) * width  + hpos;
  yval = outline{i}(:,2) * height + vpos;
  ft_plot_vector(xval, yval, 'Color','k', 'LineWidth',2)
end


% Create isolines
if strcmp(style,'iso') || strcmp(style,'surfiso')
  if ~isempty(isolines)
    contour(Xi/10,Yi,Zi,isolines,'k');
  end
end

% Plot surface
if strcmp(style,'surf') || strcmp(style,'surfiso')
  deltax = xi(2)-xi(1); % length of grid entry
  deltay = yi(2)-yi(1); % length of grid entry
  h = surf(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)), Zi, 'EdgeColor', 'none', 'FaceColor', shading);
  
  %if exist('maskimagetmp')
  %  set(h, 'facealpha', 'flat');
  %  set(h, 'alphadatamapping', 'scaled');
  %  set(h, 'alphadata', maskimagetmp);
  %end
end

% Plot filled contours
if strcmp(style,'isofill') && ~isempty(isolines)
  contourf(Xi,Yi,Zi,isolines,'k');
end



% remember the current input arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_argin     = current_argin;
previous_maskimage = maskimage;

if ~holdflag
  hold off
end

