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
% $Log: topoplotTFR.m,v $
% Revision 1.21  2009/06/17 13:44:52  roboos
% cleaned up help
%
% Revision 1.20  2008/12/16 15:31:42  sashae
% plot functions can now give cfg as output
% added checkconfig to start and end of function, configtracking possible
%
% Revision 1.19  2008/01/29 19:43:33  sashae
% added option for trial selection; plot functions now also accept data with
% repetitions (either trials or subjects), the avg is computed and plotted
% removed some old code
%
% Revision 1.18  2007/03/21 15:36:36  chrhes
% updated documentation regarding the fact that cfg.layout can also contain a
% layout structure obtained using the function prepare_layout.m
%
% Revision 1.17  2006/07/17 12:39:47  ingnie
% updated documentation
%
% Revision 1.16  2006/05/30 14:21:00  ingnie
% updated documentation
%
% Revision 1.15  2006/05/23 16:05:21  ingnie
% updated documentation
%
% Revision 1.14  2006/05/09 12:21:41  ingnie
% updated help
%
% Revision 1.13  2006/04/27 09:35:49  ingnie
% updated documentation
%
% Revision 1.12  2006/03/14 14:55:16  roboos
% changed from DOS into UNIX format
%
% Revision 1.11  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
% 

cfg=topoplotER(cfg, varargin{:});
