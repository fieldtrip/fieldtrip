function ft_plot_matrix(varargin)

% FT_PLOT_MATRIX visualizes a matrix as an image, similar to IMAGESC.
% The position, width and height can be controlled to allow multiple
% matrices (i.e. channels) to be plotted in a topographic arrangement.
%
% Use as
%   ft_plot_matrix(C, ...)
% where C is a 2 dimensional MxN matrix, or
%   ft_plot_matrix(X, Y, C, ...)
% where X and Y describe the 1xN horizontal and 1xM vertical axes
% respectively.
%
% Optional arguments should come in key-value pairs and can include
%   'clim'            = 1x2 vector with color limits (default is automatic)
%   'highlight'       = a logical matrix of size C, where 0 means that the corresponding values in C are highlighted according to the highlightstyle
%   'highlightstyle'  = can be 'saturation', 'opacity', 'outline' or 'colormix' (default = 'opacity')
%   'box'             = draw a box around the local axes, can be 'yes' or 'no'
%   'tag'             = string, the name assigned to the object. All tags with the same name can be deleted in a figure, without deleting other parts of the figure.
%
% It is possible to plot the object in a local pseudo-axis (c.f. subplot), which is specfied as follows
%   'hpos'            = horizontal position of the center of the local axes
%   'vpos'            = vertical position of the center of the local axes
%   'width'           = width of the local axes
%   'height'          = height of the local axes
%   'hlim'            = horizontal scaling limits within the local axes
%   'vlim'            = vertical scaling limits within the local axes
%
% When using a local pseudo-axis, you can plot a label next to the data
%   'label'           = string, label to be plotted at the upper left corner
%   'fontcolor'       = string, color specification (default = 'k')
%   'fontsize'        = number, sets the size of the text (default = 10)
%   'fontunits'       =
%   'fontname'        =
%   'fontweight'      =
%
% Example
%   ft_plot_matrix(randn(30,50), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0)
%
% See also FT_PLOT_VECTOR, IMAGESC, SURF

% Copyrights (C) 2009-2022, Robert Oostenveld
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
hpos            = ft_getopt(varargin, 'hpos');
vpos            = ft_getopt(varargin, 'vpos');
width           = ft_getopt(varargin, 'width');
height          = ft_getopt(varargin, 'height');
hlim            = ft_getopt(varargin, 'hlim');
vlim            = ft_getopt(varargin, 'vlim');
clim            = ft_getopt(varargin, 'clim');
highlight       = ft_getopt(varargin, 'highlight');
highlightstyle  = ft_getopt(varargin, 'highlightstyle', 'opacity');
label           = ft_getopt(varargin, 'label');
box             = ft_getopt(varargin, 'box',            false);
tag             = ft_getopt(varargin, 'tag',            '');
% these have to do with the font of the label
fontcolor       = ft_getopt(varargin, 'fontcolor', 'k'); % default is black
fontsize        = ft_getopt(varargin, 'fontsize',   get(0, 'defaulttextfontsize'));
fontname        = ft_getopt(varargin, 'fontname',   get(0, 'defaulttextfontname'));
fontweight      = ft_getopt(varargin, 'fontweight', get(0, 'defaulttextfontweight'));
fontunits       = ft_getopt(varargin, 'fontunits',  get(0, 'defaulttextfontunits'));

if ~isempty(highlight) && ~isequal(size(highlight), size(cdat))
  ft_error('the dimensions of the highlight should be identical to the dimensions of the data');
end

% axis   = ft_getopt(varargin, 'axis', false);
% style  = ft_getopt(varargin, 'style'); % FIXME

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
      ft_error('unsupported option for hlim')
  end % switch
end % if ischar

if hlim(1)==hlim(2)
  if hlim(1)==0
    % automatic scaling not possible
    hlim = [-1 1];
  else
    % adjust the scaling a bit
    hlim(1) = 0.8*hlim(1);
    hlim(2) = 1.2*hlim(2);
  end
end

if ischar(vlim)
  switch vlim
    case 'maxmin'
      vlim = [min(vdat) max(vdat)];
    case 'maxabs'
      vlim = max(abs(vdat));
      vlim = [-vlim vlim];
    otherwise
      ft_error('unsupported option for vlim')
  end % switch
end % if ischar

if vlim(1)==vlim(2)
  if vlim(1)==0
    % automatic scaling not possible
    vlim = [-1 1];
  else
    % adjust the scaling a bit
    vlim(1) = 0.8*vlim(1);
    vlim(2) = 1.2*vlim(2);
  end
end

if ischar(clim)
  switch clim
    case 'maxmin'
      clim = [min(cdat(:)) max(cdat(:))];
    case 'maxabs'
      clim = max(abs(cdat(:)));
      clim = [-clim clim];
    otherwise
      ft_error('unsupported option for clim')
  end % switch
end % if ischar

if clim(1)==clim(2)
  if clim(1)==0
    % automatic scaling not possible
    clim = [-1 1];
  else
    % adjust the scaling a bit
    clim(1) = 0.8*clim(1);
    clim(2) = 1.2*clim(2);
  end
end

% these must be floating point values and not integers, otherwise the scaling fails
hdat = double(hdat);
vdat = double(vdat);
cdat = double(cdat);
hlim = double(hlim);
vlim = double(vlim);
clim = double(clim);

if isempty(hpos)
  hpos = (hlim(1)+hlim(2))/2;
end

if isempty(vpos)
  vpos = (vlim(1)+vlim(2))/2;
end

if isempty(width)
  width = hlim(2)-hlim(1);
  if length(hdat)>1
    width = width * length(hdat)/(length(hdat)-1);
  else
    width = 1;
  end
  autowidth = true;
else
  autowidth = false;
end

if isempty(height)
  height = vlim(2)-vlim(1);
  if length(vdat)>1
    height = height * length(vdat)/(length(vdat)-1);
  else
    height = 1;
  end
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

% uimagesc is in external/fileexchange
ft_hastoolbox('fileexchange', 1);

% the uimagesc-call needs to be here to avoid calling it several times in switch-highlight
if isempty(highlight)
  h = uimagesc(hdat, vdat, cdat, clim);
  set(h,'tag',tag);
end

% the uimagesc-call needs to be inside switch-statement, otherwise 'saturation' will cause it to be called twice
if ~isempty(highlight)
  switch highlightstyle
    case 'opacity'
      % get the same scaling for 'highlight' then what we will get for cdata
      h = uimagesc(hdat, vdat, highlight);
      highlight = get(h, 'CData');
      delete(h); % this is needed because "hold on" might have been called previously, e.g. in ft_multiplotTFR
      h = uimagesc(hdat, vdat, cdat, clim);
      set(h,'tag',tag);
      if ft_platform_supports('alim')
        set(h,'AlphaData',highlight);
        set(h, 'AlphaDataMapping', 'scaled');
        alim([0 1]);
      end

    case 'saturation'
      cmap    = get(gcf, 'colormap');
      rgbcdat = cdat2rgb(cdat, cmap, clim, highlight);

      h = uimagesc(hdat, vdat, rgbcdat, clim);
      set(h,'tag',tag);

    case 'colormix'
      cmap    = get(gcf, 'colormap');
      rgbcdat = bg_rgba2rgb([1 1 1], cdat, cmap, clim, highlight, 'rampup', [0 1]);

      h = uimagesc(hdat, vdat, rgbcdat, clim);
      set(h,'tag',tag);

    case 'outline'
      % the significant voxels could be outlined with a black contour
      % plot outline
      h = uimagesc(hdat, vdat, cdat, clim);
      set(h,'tag',tag);
      [x,y] = meshgrid(hdat, vdat);
      x = interp2(x, 2); % change to 4 for round corners
      y = interp2(y, 2); % change to 4 for round corners
      contourlines = highlight==1;
      contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
      dx = mean(diff(x(1, :))); % remove for round corners
      dy = mean(diff(y(:, 1))); % remove for round corners
      holdflag = ishold;
      hold on
      contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[0 0 0],'LineWidth',2);
      if ~holdflag
        hold off % revert to the previous hold state
      end

    otherwise
      ft_error('unsupported highlightstyle')
  end % switch highlightstyle
end

if ~isempty(label)
  boxposition(1) = hpos - width/2;
  boxposition(2) = hpos + width/2;
  boxposition(3) = vpos - height/2;
  boxposition(4) = vpos + height/2;
  text(boxposition(1), boxposition(4), label, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
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
