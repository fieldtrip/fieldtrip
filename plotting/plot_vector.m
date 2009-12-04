function [varargout] = plot_vector(varargin)

% PLOT_VECTOR
%
% Use as
%   plot_vector(Y, ...)
%   plot_vector(X, Y, ...)
% where X and Y are similar as the input to the Matlab plot function.
%
% Additional options should be specified in key-value pairs and can be
%   'hpos'
%   'vpos'
%   'width'
%   'height'
%   'hlim'
%   'vlim'
%   'style'
%   'label'
%   'fontsize'
%   'axis'          can be 'yes' or 'no'
%   'box'           can be 'yes' or 'no'
%   'highlight'
%   'highlightstyle'
%   'color'
%   'linewidth'
%   'markersize'
%   'markerfacecolor'
%
%
%
%
% Example use
%   plot_vector(randn(1,100), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0)

% Copyrights (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

holdflag = ishold;
hold on

if nargin>1 && all(cellfun(@isnumeric, varargin(1:2)))
  % the function was called like plot(x, y, ...)
  hdat = varargin{1};
  vdat = varargin{2};
  varargin = varargin(3:end);
else
  % the function was called like plot(y, ...)
  vdat = varargin{1};
  if any(size(vdat)==1)
    % ensure that it is a column vector
    vdat = vdat(:);
  end
  hdat = 1:size(vdat,1);
  varargin = varargin(2:end);
end

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'hlim', 'vlim', 'style', 'label', ...
  'fontsize', 'axis', 'box','highlight','highlightstyle','color', 'linewidth','markersize','markerfacecolor'});
hpos            = keyval('hpos',               varargin);
vpos            = keyval('vpos',               varargin);
width           = keyval('width',              varargin);
height          = keyval('height',             varargin);
hlim            = keyval('hlim',               varargin); if isempty(hlim),            hlim = 'maxmin';            end
vlim            = keyval('vlim',               varargin); if isempty(vlim),            vlim = 'maxmin';            end
style           = keyval('style',              varargin); if isempty(style),           style = '-';                end
label           = keyval('label',              varargin);
fontsize        = keyval('fontsize',           varargin);
axis            = keyval('axis',               varargin); if isempty(axis),            axis = false;               end
box             = keyval('box',                varargin); if isempty(box),             box = false;                end
color           = keyval('color',              varargin);
linewidth       = keyval('linewidth',          varargin); if isempty(linewidth),       linewidth = 0.5;            end
highlight       = keyval('highlight',          varargin);
highlightstyle  = keyval('highlightstyle',     varargin); if isempty(highlightstyle),  highlightstyle = 'box';     end
markersize      = keyval('markersize',         varargin); if isempty(markersize),      markersize = 6;             end
markerfacecolor = keyval('markerfacecolor',    varargin); if isempty(markerfacecolor), markerfacecolor = 'none';   end
% convert the yes/no strings into boolean values
axis = istrue(axis);
box  = istrue(box);

% label  = keyval('label', varargin); % FIXME

if ischar(hlim)
  switch hlim
    case 'maxmin'
      hlim = [min(hdat) max(hdat)];
    case 'absmax'
      hlim = max(abs(hdat));
      hlim = [-hlim hlim];
    otherwise
      error('unsupported option for hlim')
  end % switch
end % if ischar

if ischar(vlim)
  switch vlim
    case 'maxmin'
      vlim = [min(vdat(:)) max(vdat(:))];
    case 'absmax'
      vlim = max(abs(vdat(:)));
      vlim = [-vlim vlim];
    otherwise
      error('unsupported option for vlim')
  end % switch
end % if ischar


if isempty(hpos) && ~isempty(hlim)
  hpos = (hlim(1)+hlim(2))/2;
end
if isempty(vpos) && ~isempty(vlim)
  vpos = (vlim(1)+vlim(2))/2;
end

if isempty(width) && ~isempty(hlim)
  width = hlim(2)-hlim(1);
end

if isempty(height) && ~isempty(vlim)
  height = vlim(2)-vlim(1);
end

% first shift the horizontal axis to zero
hdat = hdat - (hlim(1)+hlim(2))/2;
% then scale to length 1
if hlim(2)-hlim(1)~=0
  hdat = hdat ./ (hlim(2)-hlim(1));
else
  hdat = hdat /hlim(1);
end
% then scale to the new width
hdat = hdat .* width;
% then shift to the new horizontal position
hdat = hdat + hpos;
% first shift the vertical axis to zero
vdat = vdat - (vlim(1)+vlim(2))/2;
% then scale to length 1
if vlim(2)-vlim(1)~=0
  vdat = vdat / (vlim(2)-vlim(1));
else
  vdat = vdat /vlim(1);
end
% then scale to the new width
vdat = vdat .* height;
% then shift to the new vertical position
vdat = vdat + vpos;


% plotting lines
if isempty(color)
  h = plot(hdat, vdat, style, 'LineWidth', linewidth,'markersize',markersize,'markerfacecolor',markerfacecolor);
else
  h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color,'markersize',markersize,'markerfacecolor',markerfacecolor);
end


if ~isempty(highlight)
  switch highlightstyle
    case 'box'
      % find the sample number where the highlight begins and ends
      if ~islogical(highlight)
        highlight=logical(highlight);
        warning('converting mask to logical values')
      end
      begsample = find(diff([0 highlight 0])== 1);
      endsample = find(diff([0 highlight 0])==-1)-1;
      for i=1:length(begsample)
        begx = hdat(begsample(i));
        endx = hdat(endsample(i));
        plot_box([begx endx vpos-height/2 vpos+height/2], 'facecolor', [.6 .6 .6], 'edgecolor', 'none');
      end
      % plotting lines again, otherwise box will always be on top
      if isempty(color)
        h = plot(hdat, vdat, style, 'LineWidth', linewidth,'markersize',markersize,'markerfacecolor',markerfacecolor);
      else
        h = plot(hdat, vdat, style, 'LineWidth', linewidth, 'Color', color,'markersize',markersize,'markerfacecolor',markerfacecolor);
      end
    case 'thickness'
      % find the sample number where the highligh begins and ends
      if ~islogical(highlight)
        highlight=logical(highlight);
        warning('converting mask to logical values')
      end
      begsample = find(diff([0 highlight 0])== 1);
      endsample = find(diff([0 highlight 0])==-1)-1;
      linecolor = get(h,'Color'); % get current line color
      for i=1:length(begsample)
        hor = hdat(begsample(i):endsample(i));
        ver = vdat(begsample(i):endsample(i));
        plot(hor,ver,'linewidth',4*linewidth,'linestyle','-','Color', linecolor); % changed 3* to 4*, as 3* appeared to have no effect
      end  
    case 'saturation'
      % find the sample number where the highligh begins and ends
      if ~islogical(highlight)
        highlight=logical(highlight);
        warning('converting mask to logical values')
      end
      highlight = ~highlight; % invert the mask
      begsample = find(diff([0 highlight 0])== 1);
      endsample = find(diff([0 highlight 0])==-1)-1;
      linecolor = get(h,'Color'); % get current line color
      linecolor = (linecolor * 0.2) + 0.8; % change saturation of color
      for i=1:length(begsample)
        hor = hdat(begsample(i):endsample(i));
        ver = vdat(begsample(i):endsample(i));
        plot(hor,ver,'color',linecolor); 
      end  
    case 'opacity'
      error('unsupported highlightstyle')
    otherwise
      error('unsupported highlightstyle')
  end % switch highlightstyle
end


if ~isempty(label)
  boxposition(1) = hpos - width/2;
  boxposition(2) = hpos + width/2;
  boxposition(3) = vpos - height/2;
  boxposition(4) = vpos + height/2;
  h = text(boxposition(1), boxposition(4), label);
  if ~isempty(fontsize)
    set(h, 'Fontsize', fontsize);
  end
end

if box
  boxposition = zeros(1,4);
  % this plots a box around the original hpos/vpos with appropriate width/height
  x1 = hpos - width/2;
  x2 = hpos + width/2;
  y1 = vpos - height/2;
  y2 = vpos + height/2;

  X = [x1 x2 x2 x1 x1];
  Y = [y1 y1 y2 y2 y1];
  line(X, Y);

  %   % this plots a box around the original hpos/vpos with appropriate width/height
  %   boxposition(1) = hpos - width/2;
  %   boxposition(2) = hpos + width/2;
  %   boxposition(3) = vpos - height/2;
  %   boxposition(4) = vpos + height/2;
  %   plot_box(boxposition, 'facecolor', 'none', 'edgecolor', 'k');

  % this plots a box around the complete data
  % boxposition(1) = hlim(1);
  % boxposition(2) = hlim(2);
  % boxposition(3) = vlim(1);
  % boxposition(4) = vlim(2);
  % plot_box(boxposition, 'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim);
end

if axis
  % determine where the original [0, 0] in the data is located in the scaled and shifted axes
  x0 = interp1(hlim, hpos + [-width/2  width/2 ], 0, 'linear', 'extrap');
  y0 = interp1(vlim, vpos + [-height/2 height/2], 0, 'linear', 'extrap');

  X = [hpos-width/2  hpos+width/2];
  Y = [y0 y0];
  plot_line(X, Y);
  % str = sprintf('%g', hlim(1)); plot_text(X(1), Y(1), str);
  % str = sprintf('%g', hlim(2)); plot_text(X(2), Y(2), str);

  X = [x0 x0];
  Y = [vpos-height/2 vpos+height/2];
  plot_line(X, Y);
  % str = sprintf('%g', vlim(1)); plot_text(X(1), Y(1), str);
  % str = sprintf('%g', vlim(2)); plot_text(X(2), Y(2), str);
end

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = h;
end

if ~holdflag
  hold off
end
