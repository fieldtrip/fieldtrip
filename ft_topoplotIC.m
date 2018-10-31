function [cfg] = ft_topoplotIC(cfg, comp)

% FT_TOPOPLOTIC plots the topographic distribution of an independent
% component that was computed using the FT_COMPONENTANALYSIS function,
% as a 2-D circular view (looking down at the top of the head).
%
% Use as
%   ft_topoplotIC(cfg, comp)
% where the input comp structure should be obtained from FT_COMPONENTANALYSIS.
%
% The configuration should have the following parameters:
%   cfg.component          = field that contains the independent component(s) to be plotted as color
%   cfg.layout             = specification of the layout, see below
%
% The configuration can have the following parameters:
%   cfg.colormap           = any sized colormap, see COLORMAP
%   cfg.zlim               = plotting limits for color dimension, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.marker             = 'on', 'labels', 'numbers', 'off'
%   cfg.markersymbol       = channel marker symbol (default = 'o')
%   cfg.markercolor        = channel marker color (default = [0 0 0] (black))
%   cfg.markersize         = channel marker size (default = 2)
%   cfg.markerfontsize     = font size of channel labels (default = 8 pt)
%   cfg.highlight          = 'on', 'labels', 'numbers', 'off'
%   cfg.highlightchannel   =  Nx1 cell-array with selection of channels, or vector containing channel indices see FT_CHANNELSELECTION
%   cfg.highlightsymbol    = highlight marker symbol (default = 'o')
%   cfg.highlightcolor     = highlight marker color (default = [0 0 0] (black))
%   cfg.highlightsize      = highlight marker size (default = 6)
%   cfg.highlightfontsize  = highlight marker size (default = 8)
%   cfg.colorbar           = 'yes'
%                            'no' (default)
%                            'North'              inside plot box near top
%                            'South'              inside bottom
%                            'East'               inside right
%                            'West'               inside left
%                            'NorthOutside'       outside plot box near top
%                            'SouthOutside'       outside bottom
%                            'EastOutside'        outside right
%                            'WestOutside'        outside left
%   cfg.interplimits       = limits for interpolation (default = 'head')
%                            'electrodes' to furthest electrode
%                            'head' to edge of head
%   cfg.interpolation      = 'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
%   cfg.style              = plot style (default = 'both')
%                            'straight' colormap only
%                            'contour' contour lines only
%                            'both' (default) both colormap and contour lines
%                            'fill' constant color between lines
%                            'blank' only the head shape
%   cfg.gridscale          = scaling grid size (default = 67)
%                            determines resolution of figure
%   cfg.shading            = 'flat' 'interp' (default = 'flat')
%   cfg.comment            = string 'no' 'auto' or 'xlim' (default = 'auto')
%                            'auto': date, xparam and zparam limits are printed
%                            'xlim': only xparam limits are printed
%   cfg.commentpos         = string or two numbers, position of comment (default 'leftbottom')
%                            'lefttop' 'leftbottom' 'middletop' 'middlebottom' 'righttop' 'rightbottom'
%                            'title' to place comment as title
%                            'layout' to place comment as specified for COMNT in layout
%                            [x y] coordinates
%   cfg.title              = string or 'auto' or 'off', specify a figure
%                            title, or use 'component N' (auto) as the
%                            title
%
% The layout defines how the channels are arranged. You can specify the
% layout in a variety of ways:
%  - you can provide a pre-computed layout structure (see prepare_layout)
%  - you can give the name of an ascii layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition, i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these and the data structure contains an
% electrode or gradiometer structure, that will be used for creating a
% layout. If you want to have more fine-grained control over the layout
% of the subplots, you should create your own layout file.
%
% See also FT_COMPONENTANALYSIS, FT_REJECTCOMPONENT, FT_TOPOPLOTTFR,
% FT_SINGLEPLOTTFR, FT_MULTIPLOTTFR, FT_PREPARE_LAYOUT

% Undocumented local options:
% cfg.labeloffset (offset of labels to their marker, default = 0.005)

% Copyright (C) 2010, Donders Centre for Cognitive Neuroimaging, Arjen Stolk
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar comp
ft_preamble provenance comp
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
% this will remove all time-series information
comp = ft_checkdata(comp, 'datatype', 'comp');

% set the config defaults
cfg.title     = ft_getopt(cfg, 'title', 'auto');
cfg.parameter = ft_getopt(cfg, 'parameter', 'topo'); % needed in topoplot_common

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', 'component');
cfg = ft_checkconfig(cfg, 'allowedval', {'parameter' 'topo'});

% interactive plotting doesn't work for chan_comp dimord.
if isfield(cfg, 'interactive') && strcmp(cfg.interactive, 'yes')
  ft_warning('Interactive plotting is not supported.');
end
cfg.interactive = 'no';

% prepare the layout, this should be done only once
tmpcfg     = removefields(cfg, 'inputfile');
tmpcomp.label = comp.topolabel; % the input to ft_prepare_layout needs at least a data.label field
cfg.layout = ft_prepare_layout(tmpcfg, tmpcomp);
clear tmpcomp;

% don't show the callinfo for each separate component
tmpshowcallinfo = cfg.showcallinfo;
cfg.showcallinfo = 'no';

% create temporary variable to prevent overwriting the selected components
selcomp = cfg.component;

% make sure figure window titles are labeled appropriately, pass this onto the actual
% plotting function if we don't specify this, the window will be called
% 'ft_topoplotTFR', which is confusing to the user
cfg.funcname = mfilename;
if nargin > 1
  cfg.dataname = {inputname(2)};
  for k = 3:nargin
    cfg.dataname{end+1} = inputname(k);
  end
end

nplots = numel(selcomp);
if nplots>1
  % make multiple plots in a single figure
  nyplot = ceil(sqrt(nplots));
  nxplot = ceil(nplots./nyplot);
  for i = 1:length(selcomp)
    subplot(nxplot, nyplot, i);
    cfg.component = selcomp(i);

    % call the common function that is shared with ft_topoplotER and ft_topoplotTFR
    [cfg] = topoplot_common(cfg, comp);

    if strcmp(cfg.title, 'auto')
      title(['component ' num2str(selcomp(i))]);
    elseif ~strcmp(cfg.title, 'off')
      title(cfg.title);
    end
  end % for all components

else
  cfg.component = selcomp;

  % call the common function that is shared with ft_topoplotER and ft_topoplotTFR
  [cfg] = topoplot_common(cfg, comp);

  if strcmp(cfg.title, 'auto')
    title(['component ' num2str(selcomp)]);
  elseif ~strcmp(cfg.title, 'off')
    title(cfg.title);
  end
end

% remove this field again, it is only used for figure labels
cfg = removefields(cfg, 'funcname');

% show the callinfo for all components together
cfg.showcallinfo = tmpshowcallinfo;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous comp
ft_postamble provenance

if ~nargout
  clear cfg
end
