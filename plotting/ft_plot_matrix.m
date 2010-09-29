function ft_plot_matrix(varargin)

% FT_PLOT_MATRIX
%
% Use as
%   ft_plot_matrix(C, ...)
% where C is a 2 dimensional MxN matrix, or
%   ft_plot_matrix(X, Y, C, ...)
% where X and Y describe the 1xN horizontal and 1xM vertical axes 
% respectively.
%
% Additional options should be specified in key-value pairs and can be
%   'hpos'
%   'vpos'
%   'width'
%   'height'
%   'hlim'
%   'vlim'
%   'clim'
%   'box'                can be 'yes' or 'no'
%   'highlight' 
%   'highlightstlyle'    can be 'saturation' or 'opacity'
%   'tag'
%
% Example use
%   ft_plot_matrix(randn(30,50), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0)

% Copyrights (C) 2009, Robert Oostenveld
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

warning('on', 'MATLAB:divideByZero');

if nargin>2 && all(cellfun(@isnumeric, varargin(1:3)))
  % the function was called like imagesc(x, y, c, ...)
  hdat = varargin{1};
  vdat = varargin{2};
  cdat = varargin{3};
  varargin = varargin(4:end);
else
  % the function was called like plot(c, ...)
  cdat = varargin{1};
  vdat = 1:size(cdat,1);
  hdat = 1:size(cdat,2);
  varargin = varargin(2:end);
end

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'hlim', 'vlim', 'clim', 'box','highlight','highlightstyle','tag'});
hpos           = keyval('hpos',   varargin);
vpos           = keyval('vpos',   varargin);
width          = keyval('width',  varargin);
height         = keyval('height', varargin);
hlim           = keyval('hlim',   varargin);
vlim           = keyval('vlim',   varargin);
clim           = keyval('clim',   varargin);
box            = keyval('box',    varargin);              if isempty(box),               box = false;                    end
highlight      = keyval('highlight',       varargin);
highlightstyle = keyval('highlightstyle',  varargin);     if isempty(highlightstyle),    highlightstyle = 'opacity';     end
tag            = keyval('tag', varargin);                 if isempty(tag),               tag='';                         end

% axis   = keyval('axis',   varargin); if isempty(axis), axis = false; end
% label  = keyval('label',  varargin); % FIXME
% style  = keyval('style',  varargin); % FIXME

% convert the yes/no strings into boolean values
box  = istrue(box);

if isempty(hlim)
  hlim = 'maxmin';
end

if isempty(vlim)
  vlim = 'maxmin';
end

if isempty(clim)
  clim = 'maxmin';
end

if ischar(hlim)
  switch hlim
    case 'maxmin'
      hlim = [min(hdat) max(hdat)];
    case 'maxabs'
      hlim = max(abs(hdat));
      hlim = [-hlim hlim];
    otherwise
      error('unsupported option for hlim')
  end % switch
end % if ischar

if ischar(vlim)
  switch vlim
    case 'maxmin'
      vlim = [min(vdat) max(vdat)];
    case 'maxabs'
      vlim = max(abs(vdat));
      vlim = [-vlim vlim];
    otherwise
      error('unsupported option for vlim')
  end % switch
end % if ischar

if ischar(clim)
  switch clim
    case 'maxmin'
      clim = [min(cdat(:)) max(cdat(:))];
    case 'maxabs'
      clim = max(abs(cdat(:)));
      clim = [-clim clim];
    otherwise
      error('unsupported option for clim')
  end % switch
end % if ischar

% these must be floating point values and not integers, otherwise the scaling fails
hdat = double(hdat);
vdat = double(vdat);
cdat = double(cdat);
hlim = double(hlim);
vlim = double(vlim);
clim = double(clim);

if isempty(hpos);
  hpos = (hlim(1)+hlim(2))/2;
end

if isempty(vpos);
  vpos = (vlim(1)+vlim(2))/2;
end

if isempty(width),
  width = hlim(2)-hlim(1);
  width = width * length(hdat)/(length(hdat)-1);
  autowidth = true;
else
  autowidth = false;
end

if isempty(height),
  height = vlim(2)-vlim(1);
  height = height * length(vdat)/(length(vdat)-1);
  autoheight = true;
else
  autoheight = false;
end

% hlim
% vlim

% first shift the horizontal axis to zero
hdat = hdat - (hlim(1)+hlim(2))/2;
% then scale to length 1
hdat = hdat ./ (hlim(2)-hlim(1));
% then scale to compensate for the patch size
hdat = hdat * (length(hdat)-1)/length(hdat);
% then scale to the new width
hdat = hdat .* width;
% then shift to the new horizontal position
hdat = hdat + hpos;

% first shift the vertical axis to zero
vdat = vdat - (vlim(1)+vlim(2))/2;
% then scale to length 1
vdat = vdat ./ (vlim(2)-vlim(1));
% then scale to compensate for the patch size
vdat = vdat * (length(vdat)-1)/length(vdat);
% then scale to the new width
vdat = vdat .* height;
% then shift to the new vertical position
vdat = vdat + vpos;

% the uimagesc-call needs to be here to avoid calling it several times in switch-highlight
if isempty(highlight)
  h = uimagesc(hdat, vdat, cdat, clim);
  set(h,'tag',tag);
end

% the uimagesc-call needs to be inside switch-statement, otherwise 'saturation' will cause it to be called twice
if ~isempty(highlight)
  switch highlightstyle
    case 'opacity'
      h = uimagesc(hdat, vdat, cdat, clim);
      
      set(h,'CData',cdat); % quick fix
      
      set(h,'tag',tag);
      set(h,'AlphaData',highlight);
      set(h, 'AlphaDataMapping', 'scaled');
      alim([0 1]);
    case 'saturation'
      satmask = highlight;
      
      % Transform cdat-values to have a 0-64 range, dependent on clim
      % (think of it as the data having an exact range of min=clim(1) to max=(clim2), convert this range to 0-64)
      cdat = (cdat + -clim(1)) * (64 / (-clim(1) + clim(2))); 
      
      % Make sure NaNs are plotted as white pixels, even when using non-integer mask values
      satmask(isnan(cdat)) = 0;
      cdat(isnan(cdat)) = 32;
           
      % ind->rgb->hsv ||change saturation values||  hsv->rgb ->  plot
      rgbcdat = ind2rgb(uint8(floor(cdat)), colormap);
      hsvcdat = rgb2hsv(rgbcdat);
      hsvcdat(:,:,2) = hsvcdat(:,:,2) .* satmask;
      rgbcdatsat = hsv2rgb(hsvcdat);
      h = uimagesc(hdat, vdat, rgbcdatsat,clim);
      set(h,'tag',tag);
    case 'outline'
      % the significant voxels could be outlined with a black contour
      error('unsupported highlightstyle')
    otherwise
      error('unsupported highlightstyle')
  end % switch highlightstyle
end

if box
  boxposition = zeros(1,4);
  % this plots a box around the original hpos/vpos with appropriate width/height
  boxposition(1) = hpos - width/2;
  boxposition(2) = hpos + width/2;
  boxposition(3) = vpos - height/2;
  boxposition(4) = vpos + height/2;
  ft_plot_box(boxposition);
end
