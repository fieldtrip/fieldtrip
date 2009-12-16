function [cfg] = topoplotTFR(cfg, varargin)

% TOPOPLOTTFR plots the topographic distribution of 3-Dimensional datatypes as 
% the time-frequency representation of power or coherence that was computed
% using the FREQANALYSIS or FREQDESCRIPTIVES functions, as a 2-D circular
% view (looking down at the top of the head).
%
% Use as:
%   topoplotTFR(cfg, data)
%
% The configuration can have the following parameters:
% cfg.xparam             = first dimension in data in which a selection is made
%                         'time' or 'freq' (default depends on data.dimord)
% cfg.yparam             = second dimension in data in which a selection is made
%                         'freq' (default depends on data.dimord)
% cfg.zparam             = field that contains the data to be plotted as color 
%                         'avg', 'powspctrm' or 'cohspctrm' (default depends on data.dimord)
% cfg.xlim               = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim               = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.zlim               = 'maxmin', 'absmax' or [zmin zmax] (default = 'maxmin')
% cfg.cohrefchannel      = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline           = 'yes','no' or [time1 time2] (default = 'no'), see TIMELOCKBASELINE or FREQBASELINE
% cfg.baselinetype       = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials             = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.colormap           = any sized colormap, see COLORMAP
% cfg.marker             = 'on', 'labels', 'numbers', 'off'                    
% cfg.markersymbol       = channel marker symbol (default = 'o')
% cfg.markercolor        = channel marker color (default = [0 0 0] (black))
% cfg.markersize         = channel marker size (default = 2)
% cfg.markerfontsize     = font size of channel labels (default = 8 pt)                
% cfg.highlight          = 'on', 'labels', 'numbers', 'off'                    
% cfg.highlightchannel   =  Nx1 cell-array with selection of channels, or vector containing channel indices see CHANNELSELECTION 
% cfg.highlightsymbol    = highlight marker symbol (default = 'o')
% cfg.highlightcolor     = highlight marker color (default = [0 0 0] (black))
% cfg.highlightsize      = highlight marker size (default = 6)
% cfg.highlightfontsize  = highlight marker size (default = 8)
% cfg.colorbar           = 'yes'
%                          'no' (default)
%                          'North'              inside plot box near top
%                          'South'              inside bottom
%                          'East'               inside right
%                          'West'               inside left
%                          'NorthOutside'       outside plot box near top
%                          'SouthOutside'       outside bottom
%                          'EastOutside'        outside right
%                          'WestOutside'        outside left
% cfg.interplimits       = limits for interpolation (default = 'head')
%                          'electrodes' to furthest electrode
%                          'head' to edge of head
% cfg.interpolation      = 'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
% cfg.style              = plot style (default = 'both')
%                          'straight' colormap only
%                          'contour' contour lines only
%                          'both' (default) both colormap and contour lines
%                          'fill' constant color between lines
%                          'blank' only the head shape
% cfg.gridscale          = scaling grid size (default = 67)
%                          determines resolution of figure
% cfg.shading            = 'flat' 'interp' (default = 'flat')
% cfg.comment            = string 'no' 'auto' or 'xlim' (default = 'auto')
%                          'auto': date, xparam and zparam limits are printed
%                          'xlim': only xparam limits are printed
% cfg.commentpos         = string or two numbers, position of comment (default 'leftbottom')
%                          'lefttop' 'leftbottom' 'middletop' 'middlebottom' 'righttop' 'rightbottom'
%                          'title' to place comment as title
%                          'layout' to place comment as specified for COMNT in layout
%                          [x y] coordinates
% cfg.interactive        = Interactive plot 'yes' or 'no' (default = 'no')
%                          In a interactive plot you can select areas and produce a new
%                          interactive plot when a selected area is clicked. Multiple areas 
%                          can be selected by holding down the SHIFT key.
% cfg.layout             = specification of the layout, see below
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
% TOPOPLOTTFR calls the function TOPOPLOTER to do the plotting.
%
% See also:
%   topoplotER, singleplotTFR, multiplotTFR, prepare_layout

% Undocumented local options:
% cfg.labeloffset (offset of labels to their marker, default = 0.005)


% Copyright (C) 2005-2006, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

cfg=topoplotER(cfg, varargin{:});
