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
% cfg.xparam        = first dimension in data in which a selection is made
%                     'time' (default depends on data.dimord)
% cfg.yparam        = second dimension in data in which a selection is made
%                     'freq' (default depends on data.dimord)
% cfg.zparam        = field that contains the data to be plotted as color 
%                     'powspctrm' or 'cohspctrm' (default depends on data.dimord)
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.zlim          = 'maxmin', 'absmax' or [zmin zmax] (default = 'maxmin')
% cfg.cohrefchannel = name of reference channel for visualising coherence, can be 'gui'
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.comment       = string 'no' 'auto' or 'xlim' (default = 'auto')
%                     'auto': date, xparam, yparam and zparam limits are printed
%                     'xlim': only xparam limits are printed
% cfg.commentpos    = string or two numbers, position of comment (default 'leftbottom')
%                     'lefttop' 'leftbottom' 'middletop' 'middlebottom' 'righttop' 'rightbottom'
%                     'title' to place comment as title
%                     'layout' to place comment as specified for COMNT in layout
%                     [x y] coordinates
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas 
%                     can be selected by holding down the SHIFT key.
% cfg.layout        = specification of the layout, see below
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
% TOPOPLOTTFR calls the function TOPOPLOT to do the actual plotting. See
% the help of that function for more configuration options.
%
% See also:
%   topoplot, topoplotER, singleplotTFR, multiplotTFR, prepare_layout

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

cfg=topoplotER(cfg, varargin{:});
