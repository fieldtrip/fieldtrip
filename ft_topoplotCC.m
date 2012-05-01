function [cfg] = ft_topoplotCC(cfg, freq)

% FT_TOPOPLOTCC plots the coherence between channel pairs
%
% Use as
%  ft_topoplotCC(cfg, freq)
%
% The configuration should contain:
%   cfg.feedback    = string (default = 'textbar')
%   cfg.layout      = specification of the layout, see FT_PREPARE_LAYOUT
%   cfg.foi         = the frequency of interest which is to be plotted (default is the first frequency bin)
%   cfg.widthparam  = string, parameter to be used to control the line width
%   cfg.alphaparam  = string, parameter to be used to control the opacity
%   cfg.colorparam  = string, parameter to be used to control the line color
%
% The widthparam should be indicated in pixels, e.g. usefull numbers are 1
% and larger.
%
% The alphaparam should be indicated as opacity between 0 (fully transparent)
% and 1 (fully opaque).
%
% The default is to plot the connections as lines, but you can also use
% bidirectional arrows:
%    cfg.arrowhead    = none, stop, start, both (default = 'none')
%    cfg.arrowsize    = size of the arrow head (default = automatic)
%    cfg.arrowoffset  = amount that the arrow is shifted to the side (default = automatic)
%    cfg.arrowlength  = amount by which the length is reduced (default = 0.8)
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure. For this particular function, the input should be
% structured as a cell array.
%
% See also FT_PREPARE_LAYOUT, FT_MULTIPLOTCC, FT_CONNECTIVITYPLOT

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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar freq

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'cmbrepresentation', 'sparse');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'foi', 'layout'});

% set the defaults
if ~isfield(cfg, 'feedback'),   cfg.feedback = 'text';        end
if ~isfield(cfg, 'alphaparam'), cfg.alphaparam = [];          end
if ~isfield(cfg, 'widthparam'), cfg.widthparam = [];          end
if ~isfield(cfg, 'colorparam'), cfg.colorparam = 'cohspctrm'; end
if ~isfield(cfg, 'newfigure'),  cfg.newfigure = 'yes';        end

if ~isfield(cfg, 'arrowhead'),   cfg.arrowhead = 'none';       end % none, stop, start, both
if ~isfield(cfg, 'arrowsize'),   cfg.arrowsize = nan;          end % length of the arrow head, should be in in figure units, i.e. the same units as the layout
if ~isfield(cfg, 'arrowoffset'), cfg.arrowoffset = nan;        end % absolute, should be in in figure units, i.e. the same units as the layout
if ~isfield(cfg, 'arrowlength'), cfg.arrowlength = 0.8;        end % relative to the complete line

lay = ft_prepare_layout(cfg, freq);

beglabel = freq.labelcmb(:,1);
endlabel = freq.labelcmb(:,2);
ncmb     = size(freq.labelcmb,1);

% select the data to be used in the figure
fbin = nearest(freq.freq, cfg.foi);

if isfield(freq, cfg.widthparam)
  widthparam = freq.(cfg.widthparam)(:,fbin);
else
  widthparam = ones(ncmb,1);
end

if isfield(freq, cfg.alphaparam)
  alphaparam = freq.(cfg.alphaparam)(:,fbin);
else
  alphaparam = [];
end

if isfield(freq, cfg.colorparam)
  colorparam = freq.(cfg.colorparam)(:,:,fbin);
else
  colorparam = [];
end

if strcmp(cfg.newfigure, 'yes')
  figure
end

hold on
axis equal

% set the figure window title
funcname = mfilename();
if nargin < 2
  dataname = cfg.inputfile;
else
  dataname = inputname(2);
end
set(gcf, 'Name', sprintf('%d: %s: %s', gcf, funcname, join_str(', ',dataname)));
set(gcf, 'NumberTitle', 'off');

if isnan(cfg.arrowsize)
  % use the size of the figure to estimate a decent number
  siz = axis;
  cfg.arrowsize = (siz(2) - siz(1))/50;
  warning('using an arrowsize of %f', cfg.arrowsize);
end

if isnan(cfg.arrowoffset)
  % use the size of the figure to estimate a decent number
  siz = axis;
  cfg.arrowoffset = (siz(2) - siz(1))/100;
  warning('using an arrowoffset of %f', cfg.arrowoffset);
end

rgb  = colormap;
if ~isempty(colorparam)
  cmin = min(colorparam(:));
  cmax = max(colorparam(:));
  colorparam = (colorparam - cmin)./(cmax-cmin);
  colorparam = round(colorparam * (size(rgb,1)-1) + 1);
end

if strcmp(cfg.newfigure, 'yes')
  % also plot the position of the electrodes
  ft_plot_vector(lay.pos(:,1), lay.pos(:,2), 'style','k.');

  % also plot the outline, i.e. head shape or sulci
  if isfield(lay, 'outline')
    fprintf('solid lines indicate the outline, e.g. head shape or sulci\n');
    for i=1:length(lay.outline)
      if ~isempty(lay.outline{i})
        X = lay.outline{i}(:,1);
        Y = lay.outline{i}(:,2);
        ft_plot_line(X, Y, 'color', 'k', 'linewidth', 1.5, 'linestyle', '-');
      end
    end
  end

  % also plot the mask, i.e. global outline for masking the topoplot
  if isfield(lay, 'mask')
    fprintf('dashed lines indicate the mask for topograpic interpolation\n');
    for i=1:length(lay.mask)
      if ~isempty(lay.mask{i})
        X = lay.mask{i}(:,1);
        Y = lay.mask{i}(:,2);
        % the polygon representing the mask should be closed
        X(end+1) = X(1);
        Y(end+1) = Y(1);
        ft_plot_line(X, Y, 'color', 'k', 'linewidth', 1.5, 'linestyle', '-');
      end
    end
  end
end % if newfigure

% fix the limits for the axis
axis(axis);

ft_progress('init', cfg.feedback, 'plotting connections...');

for i=1:ncmb

  if strcmp(beglabel{i}, endlabel{i})
    % skip autocombinations
    continue
  end

  ft_progress(i/ncmb, 'plotting connection %d from %d (%s -> %s)\n', i, ncmb, beglabel{i}, endlabel{i});

  if widthparam(i)>0
    begindx = strcmp(beglabel{i}, lay.label);
    endindx = strcmp(endlabel{i}, lay.label);
    xbeg = lay.pos(begindx,1);
    ybeg = lay.pos(begindx,2);
    xend = lay.pos(endindx,1);
    yend = lay.pos(endindx,2);

    if strcmp(cfg.arrowhead, 'none')
      x = [xbeg xend]';
      y = [ybeg yend]';
      % h = line(x, y);
      h = patch(x, y, 1);
    else
      arrowbeg  = [xbeg ybeg];
      arrowend  = [xend yend];
      center    = (arrowbeg+arrowend)/2;
      direction = (arrowend - arrowbeg);
      direction = direction/norm(direction);
      offset    = [direction(2) -direction(1)];
      arrowbeg  = cfg.arrowlength * (arrowbeg-center) + center + cfg.arrowoffset * offset;
      arrowend  = cfg.arrowlength * (arrowend-center) + center + cfg.arrowoffset * offset;

      h = arrow(arrowbeg, arrowend, 'Ends', cfg.arrowhead, 'length', 0.05);

    end % if arrow

    if ~isempty(widthparam)
      set(h, 'LineWidth', widthparam(i));
    end

    if ~isempty(alphaparam)
      set(h, 'EdgeAlpha', alphaparam(i));
      set(h, 'FaceAlpha', alphaparam(i)); % for arrowheads
    end

    if ~isempty(colorparam)
      set(h, 'EdgeColor', rgb(colorparam(i),:));
      set(h, 'FaceColor', rgb(colorparam(i),:)); % for arrowheads
    end

  end
end
ft_progress('close');

% improve the fit in the axis
axis tight

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous freq


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for plotting arrows, see also fieldtrip/private/arrow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = arrow(arrowbeg, arrowend, varargin)
ends   = ft_getopt(varargin, 'ends');
length = ft_getopt(varargin, 'length'); % the length of the arrow head, in figure units
color  = [0 0 0]; % in RGB

direction = (arrowend - arrowbeg);
direction = direction/norm(direction);
offset    = [direction(2) -direction(1)];

pnt1 = arrowbeg;
pnt2 = arrowend;
h = patch([pnt1(1) pnt2(1)], [pnt1(2) pnt2(2)], color);

switch ends
  case 'stop'
    pnt1 = arrowend - length*direction + 0.4*length*offset;
    pnt2 = arrowend;
    pnt3 = arrowend - length*direction - 0.4*length*offset;
    h(end+1) = patch([pnt1(1) pnt2(1) pnt3(1)]', [pnt1(2) pnt2(2) pnt3(2)]', color);

  case 'start'
    pnt1 = arrowbeg + length*direction + 0.4*length*offset;
    pnt2 = arrowbeg;
    pnt3 = arrowbeg + length*direction - 0.4*length*offset;
    h(end+1) = patch([pnt1(1) pnt2(1) pnt3(1)]', [pnt1(2) pnt2(2) pnt3(2)]', color);

  case 'both'
    pnt1 = arrowend - length*direction + 0.4*length*offset;
    pnt2 = arrowend;
    pnt3 = arrowend - length*direction - 0.4*length*offset;
    h(end+1) = patch([pnt1(1) pnt2(1) pnt3(1)]', [pnt1(2) pnt2(2) pnt3(2)]', color);

    pnt1 = arrowbeg + length*direction + 0.4*length*offset;
    pnt2 = arrowbeg;
    pnt3 = arrowbeg + length*direction - 0.4*length*offset;
    h(end+1) = patch([pnt1(1) pnt2(1) pnt3(1)]', [pnt1(2) pnt2(2) pnt3(2)]', color);

  case 'none'
    % don't draw arrow heads
end

