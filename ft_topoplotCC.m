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
%   cfg.widthparam  = string, parameter to be used to control the line width (see below)
%   cfg.alphaparam  = string, parameter to be used to control the opacity (see below)
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
%    cfg.arrowhead    = string, 'none', 'stop', 'start', 'both' (default = 'none')
%    cfg.arrowsize    = scalar, size of the arrow head in figure units,
%                       i.e. the same units as the layout (default is automatically determined)
%    cfg.arrowoffset  = scalar, amount that the arrow is shifted to the side in figure units,
%                       i.e. the same units as the layout (default is automatically determined)
%    cfg.arrowlength  = scalar, amount by which the length is reduced relative to the complete line (default = 0.8)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure. For this particular function, the input should be
% structured as a cell-array.
%
% See also FT_PREPARE_LAYOUT, FT_MULTIPLOTCC, FT_CONNECTIVITYPLOT

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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar freq
ft_preamble provenance freq
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'cmbrepresentation', 'sparse');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'foi', 'layout'});

% set the defaults
cfg.feedback    = ft_getopt(cfg, 'feedback',    'text');
cfg.alphaparam  = ft_getopt(cfg, 'alphaparam',  []);
cfg.widthparam  = ft_getopt(cfg, 'widthparam',  []);
cfg.colorparam  = ft_getopt(cfg, 'colorparam',  'cohspctrm');
cfg.newfigure   = ft_getopt(cfg, 'newfigure',   'yes');
cfg.arrowhead   = ft_getopt(cfg, 'arrowhead',   'none');  % none, stop, start, both
cfg.arrowsize   = ft_getopt(cfg, 'arrowsize',   nan);     % length of the arrow head, should be in in figure units, i.e. the same units as the layout
cfg.arrowoffset = ft_getopt(cfg, 'arrowoffset', nan);     % absolute, should be in figure units, i.e. the same units as the layout
cfg.arrowlength = ft_getopt(cfg, 'arrowlength', 0.8);     % relative to the complete line
cfg.linestyle   = ft_getopt(cfg, 'linestyle',   []);
cfg.colormap    = ft_getopt(cfg, 'colormap',    colormap);
cfg.renderer    = ft_getopt(cfg, 'renderer',    []);      % let MATLAB decide on the default

tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'scalepos', 'elec', 'grad', 'opto', 'showcallinfo'});
lay = ft_prepare_layout(tmpcfg, freq);

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
  colorparam = freq.(cfg.colorparam)(:,fbin);
else
  colorparam = [];
end

if strcmp(cfg.newfigure, 'yes')
  figure
end

hold on
axis equal

if isnan(cfg.arrowsize)
  % use the size of the figure to estimate a decent number
  siz = axis;
  cfg.arrowsize = (siz(2) - siz(1))/50;
  ft_warning('using an arrowsize of %f', cfg.arrowsize);
end

if isnan(cfg.arrowoffset)
  % use the size of the figure to estimate a decent number
  siz = axis;
  cfg.arrowoffset = (siz(2) - siz(1))/100;
  ft_warning('using an arrowoffset of %f', cfg.arrowoffset);
end

rgb = cfg.colormap;
if ~isempty(colorparam)
  cmin = min(colorparam(:));
  cmax = max(colorparam(:));

  % this line creates a sorting vector that cause the most extreme valued arrows to be plotted last
  [srt, srtidx] = sort(abs(colorparam));

  colorparam = (colorparam - cmin)./(cmax-cmin);
  colorparam = round(colorparam * (size(rgb,1)-1) + 1);
end

if strcmp(cfg.newfigure, 'yes')
  ft_plot_layout(lay, 'label', 'no', 'box', 'off');
end

% fix the limits for the axis
axis(axis);

if ~exist('srtidx', 'var')
  srtidx = 1:ncmb;
end

ft_progress('init', cfg.feedback, 'plotting connections...');

for i=srtidx(:)'

  if strcmp(beglabel{i}, endlabel{i})
    % skip autocombinations
    continue
  end


  if widthparam(i)>0 && (isempty(alphaparam)||alphaparam(i)>0)
    ft_progress(i/ncmb, 'plotting connection %d from %d (%s -> %s)\n', i, ncmb, beglabel{i}, endlabel{i});

    begindx = strcmp(beglabel{i}, lay.label);
    endindx = strcmp(endlabel{i}, lay.label);
    xbeg = lay.pos(begindx,1);
    ybeg = lay.pos(begindx,2);
    xend = lay.pos(endindx,1);
    yend = lay.pos(endindx,2);

    if isempty(cfg.linestyle)
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
    elseif ~isempty(cfg.linestyle)

      % new style of plotting, using curved lines, this does not allow for
      % alpha mapping on the line segment
      switch cfg.linestyle
        case 'curve'
          tmp = cscvn([xbeg mean([xbeg,xend])*0.8 xend;ybeg mean([ybeg,yend])*0.8 yend]);
          pnt = fnplt(tmp);
          h   = line(pnt(1,:)', pnt(2,:)');

          if ~isempty(widthparam)
            set(h, 'LineWidth', widthparam(i));
          end
          if ~isempty(colorparam)
            set(h, 'Color', rgb(colorparam(i),:));
          end

          % deal with the arrow
          if ~strcmp(cfg.arrowhead, 'none')
            arrowbeg  = pnt(:,1)';
            arrowend  = pnt(:,end)';
            directionbeg = (arrowbeg - pnt(:,2)');
            directionend = (arrowend - pnt(:,end-1)');

            directionbeg = directionbeg/norm(directionbeg);
            directionend = directionend/norm(directionend);

            switch cfg.arrowhead
              case 'stop'
                pnt1 = arrowend - 0.05*directionend + 0.02*[directionend(2) -directionend(1)];
                pnt2 = arrowend;
                pnt3 = arrowend - 0.05*directionend - 0.02*[directionend(2) -directionend(1)];
                h_arrow = patch([pnt1(1) pnt2(1) pnt3(1)]', [pnt1(2) pnt2(2) pnt3(2)]', [0 0 0]);

              case 'start'
                pnt1 = arrowbeg - 0.05*directionbeg + 0.02*[directionbeg(2) -directionbeg(1)];
                pnt2 = arrowbeg;
                pnt3 = arrowbeg - 0.05*directionbeg - 0.02*[directionbeg(2) -directionbeg(1)];
                h_arrow = patch([pnt1(1) pnt2(1) pnt3(1)]', [pnt1(2) pnt2(2) pnt3(2)]', [0 0 0]);

              case 'both'
                pnt1 = arrowbeg - 0.05*directionbeg + 0.02*[directionbeg(2) -directionbeg(1)];
                pnt2 = arrowbeg;
                pnt3 = arrowbeg - 0.05*directionbeg - 0.02*[directionbeg(2) -directionbeg(1)];
                h_arrow(1) = patch([pnt1(1) pnt2(1) pnt3(1)]', [pnt1(2) pnt2(2) pnt3(2)]', [0 0 0]);

                pnt1 = arrowend - 0.05*directionend + 0.02*[directionend(2) -directionend(1)];
                pnt2 = arrowend;
                pnt3 = arrowend - 0.05*directionend - 0.02*[directionend(2) -directionend(1)];
                h_arrow(2) = patch([pnt1(1) pnt2(1) pnt3(1)]', [pnt1(2) pnt2(2) pnt3(2)]', [0 0 0]);

            end
          else
            h_arrow = [];
          end

          if ~isempty(widthparam)
            set(h, 'LineWidth', widthparam(i));
            if ~isempty(h_arrow)
              set(h_arrow, 'LineWidth', widthparam(i));
            end
          end
          if ~isempty(colorparam)
            set(h, 'Color', rgb(colorparam(i),:));
            if ~isempty(h_arrow)
              set(h_arrow, 'Edgecolor', rgb(colorparam(i),:));
              set(h_arrow, 'Facecolor', rgb(colorparam(i),:));
            end
          end

          if ~isempty(alphaparam)
            if ~isempty(h_arrow)
              set(h_arrow, 'EdgeAlpha', alphaparam(i));
              set(h_arrow, 'FaceAlpha', alphaparam(i)); % for arrowheads
            end
          end


        otherwise
          ft_error('unsupported linestyle specified');
      end

    end % if empty linestyle
  end 
end % for srtindx

ft_progress('close');

% this is needed for the figure title
if isfield(cfg, 'dataname') && ~isempty(cfg.dataname)
  dataname = cfg.dataname;
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  dataname = cfg.inputfile;
elseif nargin>1
  dataname = arrayfun(@inputname, 2:nargin, 'UniformOutput', false);
else
  dataname = {};
end

% set the figure window title
if ~isempty(dataname)
  set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), mfilename, join_str(', ', dataname)));
else
  set(gcf, 'Name', sprintf('%d: %s', double(gcf), mfilename));
end
set(gcf, 'NumberTitle', 'off');

% set renderer if specified
if ~isempty(cfg.renderer)
  set(gcf, 'renderer', cfg.renderer)
end

% improve the fit in the axis
axis tight

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous freq
ft_postamble provenance
ft_postamble savefig

% add a menu to the figure, but only if the current figure does not have subplots
menu_fieldtrip(gcf, cfg, false);

if ~ft_nargout
  % don't return anything
  clear cfg
end

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
