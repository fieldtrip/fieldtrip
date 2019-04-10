function varargout=spm(varargin)
% SPM: Statistical Parametric Mapping (startup function)
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  
% \__ \ )___/ )    (   Statistical Parametric Mapping
% (___/(__)  (_/\/\_)  SPM - http://www.fil.ion.ucl.ac.uk/spm
%_______________________________________________________________________
%
% SPM (Statistical Parametric Mapping) is a package for the analysis
% functional brain mapping experiments. It is the in-house package of
% the Wellcome Department of Cognitive Neurology, and is available to
% the scientific community as copyright freeware under the terms of the
% GNU General Public Licence.
% 
% Theoretical, computational and other details of the package are
% available in SPM's "Help" facility. This can be launched from the
% main SPM Menu window using the "Help" button, or directly from the
% command line using the command `spm_help`.
%
% Details of this release are availiable via the "About SPM" help topic
% (file spm.man), accessible from the SPM splash screen.  (Or type
% `spm_help spm.man` in the MatLab command window)
% 
% This spm function initialises the default parameters, and displays a
% splash screen with buttons leading to the PET(SPECT) & fMRI
% modalities Alternatively, `spm('pet')` and `spm('fmri')`
% (equivalently `spm pet` and `spm mri`) lead directly to the respective
% modality interfaces.
%
% Once the modality is chosen, (and it can be toggled mid-session) the
% SPM user interface is displayed. This provides a constant visual
% environment in which data analysis is implemented. The layout has
% been designed to be simple and at the same time show all the
% facilities that are available. The interface consists of three
% windows: A menu window with pushbuttons for the SPM routines (each
% button has a 'CallBack' string which launches the appropriate
% function/script); A blank panel used for interaction with the user;
% And a graphics figure with various editing and print facilities (see
% spm_figure.m). (These windows are 'Tag'ged 'Menu', 'Interactive', and
% 'Graphics' respectively, and should be referred to by their tags
% rather than their figure numbers.)
%
% Further interaction with the user is (mainly) via questioning in the
% 'Interactive' window (managed by spm_input), and file selection
% (managed by spm_get). See the help on spm_input.m and spm_get.m for
% details on using these functions.
%
% If a "message of the day" file named spm_motd.man exists in the SPM
% directory (alongside spm.m) then it is displayed in the Graphics
% window on startup.
%
% Arguments to this routine (spm.m) lead to various setup facilities,
% mainly of use to SPM power users and programmers. See programmers
% FORMAT & help in the main body of spm.m
%
%_______________________________________________________________________
% SPM is developed by members and collaborators of the
%                             Wellcome Department of Cognitive Neurology

%-SCCS ID and authorship of this program...
%-----------------------------------------------------------------------
% @(#)spm.m	2.87 Andrew Holmes 03/02/28

%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take. Recall   )
%( MatLab's command-function duality: `spm Welcome` is equivalent to   )
%( `spm('Welcome')`.                                                   )
%
% FORMAT spm
% Defaults to spm('Welcome')
%
% FORMAT spm('Welcome')
% Clears command window, deletes all figures, prints welcome banner and
% splash screen, sets window defaults.
%
% FORMAT spm('AsciiWelcome')
% Prints ASCII welcome banner in MatLab command window.
%
% FORMAT spm('PET') spm('FMRI')
% Closes all windows and draws new Menu, Interactive, and Graphics
% windows for an SPM session. The buttons in the Menu window launch the
% main analysis routines.
%
% FORMAT Fmenu = spm('CreateMenuWin',Vis)
% Creates SPM menu window, 'Tag'ged 'Menu'
% F   - handle of figure created
% Vis - Visibility, 'on' or 'off'
%
% Finter = FORMAT spm('CreateIntWin',Vis)
% Creates an SPM Interactive window, 'Tag'ged 'Interactive'
% F   - handle of figure created
% Vis - Visibility, 'on' or 'off'
%
% FORMAT spm('ChMod',Modality)
% Changes modality of SPM: Currently SPM supports PET & MRI modalities,
% each of which have a slightly different Menu window and different
% defaults. This function switches to the specified modality, setting
% defaults and displaying the relevant buttons.
%
% FORMAT spm('defaults',Modality)
% Sets default global variables for the specified modality.
%
% FORMAT [Modality,ModNum]=spm('CheckModality',Modality)
% Checks the specified modality against those supported, returns
% upper(Modality) and the Modality number, it's position in the list of
% supported Modalities.
%
% FORMAT WS=spm('WinScale')
% Returns ratios of current display dimensions to that of a 1152 x 900
% Sun display. WS=[Xratio,Yratio,Xratio,Yratio]. Used for scaling other
% GUI elements.
% (Function duplicated in spm_figure.m, repeated to reduce inter-dependencies.)
%
% FORMAT [FS,sf] = spm('FontSize',FS)
% FORMAT [FS,sf] = spm('FontSizes',FS)
% Returns fontsizes FS scaled for the current display.
% FORMAT sf = spm('FontScale')
% Returns font scaling factor
% FS     - (vector of) Font sizes to scale [default [1:36]]
% sf     - font scaling factor (FS(out) = floor(FS(in)*sf)
%
% Rect = spm('WinSize',Win,raw)
% Returns sizes and positions for SPM windows.
% Win  - 'Menu', 'Interactive', 'Graphics', or '0'
%      -  Window whose position is required. Only first character is
%         examined. '0' returns size of root workspace.
% raw  - If specified, then positions are for a 1152 x 900 Sun display.
%        Otherwise the positions are scaled for the current display.
%
% FORMAT SPMdir=spm('Dir',Mfile)
% Returns the directory containing the version of spm in use,
% identified as the first in MATLABPATH containing the Mfile spm (this
% file) (or Mfile if specified).
%
% FORMAT [v,c]=spm('Ver',Mfile,ReDo,Cache,Con)
% Returns the current version (v) & copyright notice, extracted from
% the top line of the Contents.m file in the directory containing the
% currently used file Mfile (defaults on omission or empty to 'spm').
%
%-The version and copyright information are saved in a global
% variable called [upper(spm_str_manip(Mfile,'rt')),'_VER'], as a
% structure with fields 'v' and 'c'. This enables repeat use without
% recomputation.
%
%-If Con [default (missing or empty) 1] is false, then the version
% information is extracted from Mfile itself, rather than the
% Contents.m file in the same directory. When using a Contents.m file,
% the first line is read. For other files, the second line (the H1 help
% line) is used. This is for consistency with MatLab's ver and help
% commands respectively. (This functionality enables toolboxes to be
% identified by a function rather than a Contents.m file, allowing
% installation in a directory which already has a Contents.m file.)
%
%-If Cache [default (missing or empty) 1] is true, then the version and
% copyright information cached in the global variable
% [upper(Mfile),'_VER'], as a structure with fields 'v' and 'c'. This
% enables repeat use without recomputation.
%
%-If ReDo [default (missing or empty) 0] is true, then the version and
% copyright information are recomputed (regardless of any stored global
% data).
%
% FORMAT xTB = spm('TBs')
% Identifies installed SPM toolboxes: SPM toolboxes are defined as the
% contents of sub-directories of fullfile(spm('Dir'),'toolbox') - the
% SPM toolbox installation directory. For SPM to pick a toolbox up,
% there must be a single mfile in the directory whose name ends with
% the toolbox directory name. (I.e. A toolbox called "test" would be in
% the "test" subdirectory of spm('Dir'), with a single file named
% *test.m.) This M-file is regarded as the launch file for the
% toolbox.
% xTB      - structure array containing toolbox definitions
% xTB.name - name of toolbox (taken as toolbox directory name)
% xTB.prog - launch program for toolbox
% xTB.dir  - toolbox directory
%
% FORMAT spm('TBlaunch',xTB,i)
% Launch a toolbox, prepending TBdir to path if necessary
% xTB      - toolbox definition structure (i.e. from spm('TBs')
% xTB.name - name of toolbox
% xTB.prog - name of program to launch toolbox
% xTB.dir  - toolbox directory (prepended to path if not on path)
%
% FORMAT [c,cName] = spm('Colour')
% Returns the RGB triple and a description for the current en-vogue SPM
% colour, the background colour for the Menu and Help windows.
%
% FORMAT [v1,v2,...] = spm('GetGlobal',name1,name2,...)
% Returns values of global variables (without declaring them global)
% name1, name2,... - name strings of desired globals
% a1, a2,...       - corresponding values of global variables with given names
%                    ([] is returned as value if global variable doesn't exist)
%
% FORMAT CmdLine = spm('CmdLine',CmdLine)
% Command line SPM usage?
% CmdLine (input)  - CmdLine preference
%                    [defaults (missing or empty) to global defaults.cmdline,]
%                    [if it exists, or 0 (GUI) otherwise.                    ]
% CmdLine (output) - true if global CmdLine if true,
%                    or if on a terminal with no support for graphics windows.
%
% FORMAT v = spm('MLver')
% Returns MatLab version, truncated to major & minor revision numbers
%
% FORMAT spm('SetCmdWinLabel',WinStripe,IconLabel)
% Sets the names on the headers and icons of Sun command tools.
% WinStripe defaults to a summary line identifying the user, host and
% MatLab version; IconLabel to 'MatLab'.
%
% FORMAT spm('PopUpCB',h)
% Callback handler for PopUp UI menus with multiple callbacks as cellstr UserData
%
% FORMAT str = spm('GetUser',fmt)
% Returns current users login name, extracted from the hosting environment
% fmt   - format string: If USER is defined then sprintf(fmt,USER) is returned
%
% FORMAT spm('Beep')
% plays the keyboard beep!
%
% FORMAT spm('time')
% Returns the current time and date as hh:mm dd/mm/yyyy
%
% FORMAT spm('Pointer',Pointer)
% Changes pointer on all SPM (HandleVisible) windows to type Pointer
% Pointer defaults to 'Arrow'. Robust to absence of windows
%
% FORMAT h = spm('alert',Message,Title,CmdLine,wait)
% FORMAT h = spm('alert"',Message,Title,CmdLine,wait)
% FORMAT h = spm('alert*',Message,Title,CmdLine,wait)
% FORMAT h = spm('alert!',Message,Title,CmdLine,wait)
% Displays an alert, either in a GUI msgbox, or as text in the command window.
%  ( 'alert"' uses the 'help' msgbox icon, 'alert*' the )
%  ( 'error' icon, 'alert!' the 'warn' icon             )
% Message - string (or cellstr) containing message to print
% Title   - title string for alert
% CmdLine - CmdLine preference [default spm('CmdLine')]
%         - If CmdLine is complex, then a CmdLine alert is always used,
%           possibly in addition to a msgbox (the latter according
%           to spm('CmdLine').)
% wait    - if true, waits until user dismisses GUI / confirms text alert
%           [default 0] (if doing both GUI & text, waits on GUI alert)
% h       - handle of msgbox created, empty if CmdLine used
%
% FORMAT SPMid = spm('FnBanner', Fn,FnV)
% Prints a function start banner, for version FnV of function Fn, & datestamps
% FORMAT SPMid = spm('SFnBanner',Fn,FnV)
% Prints a sub-function start banner
% FORMAT SPMid = spm('SSFnBanner',Fn,FnV)
% Prints a sub-sub-function start banner
% Fn    - Function name (string)
% FnV   - Function version (string)
% SPMid - ID string: [SPMver: Fn (FnV)] 
%
% FORMAT [Finter,Fgraph,CmdLine] = spm('FnUIsetup',Iname,bGX,CmdLine)
% Robust UIsetup procedure for functions:
%   Returns handles of 'Interactive' and 'Graphics' figures.
%   Creates 'Interactive' figure if ~CmdLine, creates 'Graphics' figure if bGX.
% Iname   - Name for 'Interactive' window
% bGX     - Need a Graphics window? [default 1]
% CmdLine - CommandLine usage? [default spm('CmdLine')]
% Finter  - handle of 'Interactive' figure
% Fgraph  - handle of 'Graphics' figure
% CmdLine - CommandLine usage?
%
% FORMAT F = spm('FigName',Iname,F,CmdLine)
% Set name of figure F to "SPMver (User): Iname" if ~CmdLine
% Robust to absence of figure.
% Iname      - Name for figure
% F (input)  - Handle (or 'Tag') of figure to name [default 'Interactive']
% CmdLine    - CommandLine usage? [default spm('CmdLine')]
% F (output) - Handle of figure named
%
% FORMAT spm('GUI_FileDelete')
% CallBack for GUI for file deletion, using spm_get and confirmation dialogs
%
% FORMAT Fs = spm('Show')
% Opens all SPM figure windows (with HandleVisibility) using `figure`.
%   Maintains current figure.
% Fs - vector containing all HandleVisible figures (i.e. get(0,'Children'))
%
% FORMAT spm('Clear',Finter, Fgraph)
% Clears and resets SPM-GUI, clears and timestamps MatLab command window
% Finter  - handle or 'Tag' of 'Interactive' figure [default 'Interactive']
% Fgraph  - handle or 'Tag' of 'Graphics' figure [default 'Graphics']
%
% FORMAT spm('Help',varargin)
% Merely a gateway to spm_help(varargin) - so you can type "spm help"
% 
%_______________________________________________________________________


%-Parameters
%-----------------------------------------------------------------------
Modalities = str2mat('PET','FMRI');

%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action='Welcome'; else, Action = varargin{1}; end


%=======================================================================
switch lower(Action), case 'welcome'             %-Welcome splash screen
%=======================================================================

spm_defaults;
if isfield(defaults,'modality'), spm(defaults.modality); return; end;

%-Open startup window, set window defaults
%-----------------------------------------------------------------------
S  = get(0,'ScreenSize');
if all(S==1), error('Can''t open any graphics windows...'), end
[SPMver,SPMc] = spm('Ver','',1);
PF = spm_platform('fonts');

F = figure('IntegerHandle','off',...
	'Name',sprintf('%s%s',spm('ver'),spm('GetUser',' (%s)')),...
	'NumberTitle','off',...
	'Tag','Welcome',...
	'Position',[S(3)/2-300,S(4)/2-140,500,280],...
	'Resize','off',...
	'Pointer','Watch',...
	'Color',[1 1 1]*.8,...
	'MenuBar','none',...
	'DefaultUicontrolFontName',PF.helvetica,...
	'HandleVisibility','off',...
	'Visible','off');

%-Text
%-----------------------------------------------------------------------
hA = axes('Parent',F,'Position',[0 0 100/500 280/280],'Visible','Off');
text(0.5,0.5,'SPM',...
	'Parent',hA,...
	'FontName',PF.times,'FontSize',96,...
	'FontAngle','Italic','FontWeight','Bold',...
	'Rotation',90,...
	'VerticalAlignment','Middle','HorizontalAlignment','Center',...
	'Color',[1 1 1]*.6);

uicontrol(F,'Style','Text','Position',[110 245 390 030],...
	'String','Statistical Parametric Mapping',...
	'FontName',PF.times,'FontSize',18,'FontAngle','Italic',...
	'FontWeight','Bold',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8);


uicontrol(F,'Style','Frame','Position',[110 130 380 115]);
uicontrol(F,'Style','Frame','Position',[110 015 380 087]);

uicontrol(F,'Style','Text','String',SPMver,...
	'ToolTipString','by the FIL methods group',...
	'Position',[112 200 376 030],...
	'FontName',PF.times,'FontSize',18,'FontWeight','Bold',...
	'ForegroundColor','b')

uicontrol(F,'Style','Text','Position',[112 175 376 020],...
	'String','developed by members and collaborators of',...
	'ToolTipString','',...
	'FontName',PF.times,'FontSize',10,'FontAngle','Italic')

uicontrol(F,'Style','Text','Position',[112 160 376 020],...
	'String','The Wellcome Department of Imaging Neuroscience',...
	'ToolTipString','',...
	'FontName',PF.times,'FontSize',12,'FontWeight','Bold')

uicontrol(F,'Style','Text', 'Position',[112 140 376 020],...
	'String','Institute of Neurology, University College London',...
	'ToolTipString','',...
	'FontName',PF.times,'FontSize',12)

uicontrol(F,'Style','Text','String',SPMc,'ToolTipString',SPMc,...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontName',PF.times,'FontSize',10,...
	'HorizontalAlignment','center',...
	'Position',[110 003 380 010])

%-Objects with Callbacks - PET, fMRI, About SPM, SPMweb
%-----------------------------------------------------------------------
set(F,'DefaultUicontrolFontSize',12,'DefaultUicontrolInterruptible','on')
uicontrol(F,'String','PET and SPECT',...
	'ToolTipString',...
	'launch SPM-GUI in PET/SPECT modality (or type "spm pet" in MatLab)',...
	'Position',[140 061 150 030],...
	'CallBack','delete(gcbf),clear all,spm(''PET'')',...
	'ForegroundColor',[0 1 1])

uicontrol(F,'String','fMRI time-series',...
	'ToolTipString',...
	'launch SPM-GUI in fMRI modality (or type "spm fmri" in MatLab)',...
	'Position',[310 061 150 030],...
	'CallBack','delete(gcbf),clear all,spm(''FMRI'')',...
	'ForegroundColor',[0 1 1])

uicontrol(F,'String','About SPM',...
	'ToolTipString','launch SPMhelp browser - about SPM',...
	'Position',[140 025 100 030],...
	'CallBack','spm_help(''spm.man'')',...
	'ForegroundColor','g')

uicontrol(F,'String','SPMweb',...
	'FontWeight','Bold','FontName',PF.courier,...
	'ToolTipString',...
	'launch web browser - http://www.fil.ion.ucl.ac.uk/spm',...
	'Position',[250 025 100 030],...
	'CallBack',['set(gcbf,''Pointer'',''Watch''),',...
			'web(''http://www.fil.ion.ucl.ac.uk/spm'');',...
			'set(gcbf,''Pointer'',''Arrow'')'],...
	'ForegroundColor','k')

uicontrol(F,'String','Quit',...
	'ToolTipString','close this splash screen',...
	'Position',[360 025 100 030],...
	'CallBack','delete(gcbf)',...
	'Interruptible','off',...
	'ForegroundColor','r')

set(F,'Pointer','Arrow','Visible','on')


%=======================================================================
case 'asciiwelcome'                           %-ASCII SPM banner welcome
%=======================================================================
disp( ' ___  ____  __  __                                                  ')
disp( '/ __)(  _ \(  \/  )                                                 ')
disp( '\__ \ )___/ )    (   Statistical Parametric Mapping                 ')
disp(['(___/(__)  (_/\/\_)  ',spm('Ver'),' - http://www.fil.ion.ucl.ac.uk/spm'])
fprintf('\n')


%=======================================================================
case {'pet','fmri'}             %-Initialise SPM in PET or fMRI modality
%=======================================================================
% spm(Modality)

%-Turn on warning messages for debugging - disabled as syntax for Matlab 6.5 has changed
%warning always, warning backtrace

%-Initialisation and workspace canonicalisation
%-----------------------------------------------------------------------
clc, spm('SetCmdWinLabel')
spm('AsciiWelcome'),			fprintf('\n\nInitialising SPM')
Modality = upper(Action);					fprintf('.')
delete(get(0,'Children')),					fprintf('.')

%-Draw SPM windows
%-----------------------------------------------------------------------
Fmenu  = spm('CreateMenuWin','off');				fprintf('.')
Finter = spm('CreateIntWin','off');				fprintf('.')
spm_figure('WaterMark',Finter,spm('Ver'),'',45),		fprintf('.')
Fgraph = spm_figure('Create','Graphics','Graphics','off');	fprintf('.')

Fmotd  = fullfile(spm('Dir'),'spm_motd.man');
if exist(Fmotd), spm_help('!Disp',Fmotd,'',Fgraph,spm('Ver')); end
								fprintf('.')

%-Load startup global defaults
%-----------------------------------------------------------------------
spm_defaults,							fprintf('.')

%-Setup for current modality
%-----------------------------------------------------------------------
spm('ChMod',Modality),						fprintf('.')

%-Reveal windows
%-----------------------------------------------------------------------
set([Fmenu,Finter,Fgraph],'Visible','on')
							fprintf('done\n\n')

%-Print present working directory
%-----------------------------------------------------------------------
fprintf('SPM present working directory:\n\t%s\n',pwd)


%=======================================================================
case 'createmenuwin'                              %-Draw SPM menu window
%=======================================================================
% Fmenu = spm('CreateMenuWin',Vis)
if nargin<2, Vis='on'; else, Vis=varargin{2}; end

%-Close any existing 'Menu' 'Tag'ged windows
delete(spm_figure('FindWin','Menu'))

%-Get size and scalings and create Menu window
%-----------------------------------------------------------------------
WS   = spm('WinScale');				%-Window scaling factors
FS   = spm('FontSizes');			%-Scaled font sizes
PF   = spm_platform('fonts');			%-Font names (for this platform)
Rect = spm('WinSize','Menu','raw').*WS;		%-Menu window rectangle

[SPMver,SPMc] = spm('Ver','',1,1);

Fmenu = figure('IntegerHandle','off',...
	'Name',sprintf('%s%s',spm('ver'),spm('GetUser',' (%s)')),...
	'NumberTitle','off',...
	'Tag','Menu',...
	'Position',Rect,...
	'Resize','off',...
	'Color',[1 1 1]*.8,...
	'UserData',struct('SPMver',SPMver,'SPMc',SPMc),...
	'MenuBar','none',...
	'DefaultTextFontName',PF.helvetica,...
	'DefaultTextFontSize',FS(10),...
	'DefaultUicontrolFontName',PF.helvetica,...
	'DefaultUicontrolFontSize',FS(12),...
	'DefaultUicontrolInterruptible','on',...
	'Renderer','zbuffer',...
	'Visible','off');


%-Frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 145 380 295].*WS)

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 340 360 90].*WS)
uicontrol(Fmenu,'Style','Text',...
	'String','Spatial pre-processing',...
	'Position',[025 410 350 020].*WS,...
	'ForegroundColor','w',...
	'FontAngle','Italic',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 245 360 90].*WS)
uicontrol(Fmenu,'Style','Text',...
	'String','Model specification & parameter estimation',...
	'Position',[025 315 350 020].*WS,...
	'ForegroundColor','w',...
	'FontAngle','Italic',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 200 360 40].*WS)
uicontrol(Fmenu,'Style','Text',...
	'String','Inference',...
	'Position',[025 220 350 020].*WS,...
	'ForegroundColor','w',...
	'FontAngle','Italic',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 155 360 40].*WS)

uicontrol(Fmenu,'Style','Text',...
	'String','SPM for PET/SPECT',...
	'ToolTipString','modality & defaults set for PET/SPECT',...
	'ForegroundColor',[1 1 1]*.6,...
	'BackgroundColor',[1 1 1]*.8,...
	'FontAngle','Italic',...
	'HorizontalAlignment','center',...
	'Position',[020 122 360 020].*WS,...
	'Tag','PET',...
	'Visible','off')
uicontrol(Fmenu,'Style','Text',...
	'String','SPM for functional MRI',...
	'ToolTipString','modality & defaults set for fMRI',...
	'ForegroundColor',[1 1 1]*.6,...
	'BackgroundColor',[1 1 1]*.8,...
	'FontAngle','Italic',...
	'HorizontalAlignment','center',...
	'Position',[020 122 360 020].*WS,...
	'Tag','FMRI',...
	'Visible','off')


% Utilities frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame',...
	'BackgroundColor',spm('Colour'),...
	'Position',[010 010 380 112].*WS)

uicontrol(Fmenu,'Style','Text','String',SPMc,...
	'ToolTipString',SPMc,...
	'ForegroundColor',[1 1 1]*.6,...
	'BackgroundColor',[1 1 1]*.8,...
	'FontName',PF.times,...
	'FontSize',FS(10),...
	'HorizontalAlignment','center',...
	'Position',[020 002 360 008].*WS)

%-Objects with Callbacks - main spm_*_ui.m routines
%=======================================================================

%-Spatial
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','popup',...
	'String','Realign|Realign & Unwarp',...
	'Position',[035 380 100 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','realignment',...
	'UserData','spm_realign_ui',...
	'CallBack','spm_realign_ui;',...
	'Tag','FMRI','Visible','off');                                                                              

uicontrol(Fmenu,'String','Realign',...
	'Position',[040 380 080 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','realignment',...
	'UserData','spm_realign_ui',...
	'CallBack','spm_realign_ui;',...
	'Tag','PET','Visible','off');

uicontrol(Fmenu,'String','Normalize',...
	'Position',[150 370 100 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','spatial normalization',...
	'CallBack','spm_normalise_ui;',...
	'Tag','PET',...
	'UserData','spm_normalise_ui',...
	'Visible','off');

uicontrol(Fmenu,'String','Normalize',...
	'Position',[150 345 100 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','spatial normalisation',...
	'CallBack','spm_normalise_ui;',...
	'UserData','spm_normalise_ui',...
	'Tag','FMRI', 'Visible','off');

uicontrol(Fmenu,'String','Slice timing',...
	'Position',[150 380 100 030].*WS,...
	'FontSize',FS(10),'ToolTipString',...
	'correct slice acquisition times',...
	'CallBack','spm_slice_timing;',...
	'UserData','spm_slice_timing',...
	'Tag','FMRI',...
	'Visible','off');;

uicontrol(Fmenu,'String','Smooth',...
	'Position',[280 380 080 030].*WS,...
	'ToolTipString','spatial smoothing with Gaussian kernel',...
	'FontSize',FS(10),...
	'UserData','spm_smooth_ui',...
	'CallBack','spm_smooth_ui;');

uicontrol(Fmenu,'String','Coregister',...
	'Position',[040 345 080 030].*WS,...
	'ToolTipString','co-register images from disparate modalities',...
	'FontSize',FS(10),...
	'UserData','spm_coreg_ui',...
	'CallBack','spm_coreg_ui;');

uicontrol(Fmenu,'String','Segment',...
	'Position',[280 345 080 030].*WS,...
	'ToolTipString','segment',...
	'FontSize',FS(10),...
	'UserData','spm_segment_ui',...
	'CallBack','spm_segment_ui;');

%-Statistical
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Basic models',...
	'Position',[035 285 110 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','basic statistical models',...
	'UserData','spm_spm_ui',...
	'CallBack','SPM = spm_spm_ui(''cfg'',spm_spm_ui(''DesDefs_Stats''));');

uicontrol(Fmenu,'String','PET',...
	'Position',[150 285 95 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','GLM setup for PET',...
	'CallBack','SPM = spm_spm_ui(''cfg'',spm_spm_ui(''DesDefs_PET''));',...
	'Visible','off',...
	'Tag','PET',...
	'UserData','spm_spm_ui',...
	'Enable','on');

uicontrol(Fmenu,'String','fMRI',...
	'Position',[150 285 95 030].*WS,...
	'ToolTipString',['GLM setup for fMRI time series'],...
	'FontSize',FS(10),...
	'CallBack','SPM = spm_fmri_spm_ui;',...
	'Visible','off',...
	'UserData','spm_fmri_spm_ui',...
	'Tag','FMRI');

uicontrol(Fmenu,'String','Review design',...
	'Position',[250 285 110 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','review a specified model',...
	'UserData','spm_DesRep',...
	'CallBack','spm pointer watch, spm_DesRep; spm pointer arrow');


uicontrol(Fmenu,'String','Estimate',...
	'Position',[035 250 160 030].*WS,...
	'ToolTipString','estimate (WLS) the parameters of a specified model',...
	'FontSize',FS(10),...
	'UserData','spm_spm',...
	'CallBack','SPM = spm_spm;');


uicontrol(Fmenu,'String','-> Bayesian',...
	'Position',[200 250 160 030].*WS,...
	'ToolTipString','Bayesian or conditional estimators',...
	'FontSize',FS(10),...
	'UserData','spm_spm_Bayes',...
	'CallBack','spm_spm_Bayes',...
	'ForeGroundColor',[1 1 1]/2);

uicontrol(Fmenu,'String','Results',...
	'Position',[130 205 130 030].*WS,...
	'ToolTipString','Inference and regional responses etc.',...
	'FontSize',FS(10),...
	'UserData','spm_results_ui',...
	'CallBack','[hReg,xSPM,SPM] = spm_results_ui;');

uicontrol(Fmenu,'String','Dynamic Causal Modelling',...
	'Position',[60 160 280 030].*WS,...
	'ToolTipString','Specify and estimate a dynamic causal model',...
	'FontSize',FS(10),...
	'Visible','off',...
	'Tag','FMRI',...
	'UserData','spm_dcm_ui',...
	'CallBack','spm_dcm_ui',...
	'ForeGroundColor',[1 1 1]/2);


%-Utility buttons (first line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Display',...
	'Position',[020 088 082 024].*WS,...
	'ToolTipString','orthogonal sections',...
	'FontSize',FS(9),...
	'UserData','spm_image',...
	'CallBack','spm_image')

uicontrol(Fmenu,'String','Check Reg',...
	'Position',[112 088 083 024].*WS,...
	'ToolTipString','check image registration',...
	'FontSize',FS(9),...
	'UserData','spm_check_registration',...
	'CallBack','spm_check_registration;')

uicontrol(Fmenu,'Style','PopUp',...
	'String','Render...|Display|Xtract Brain',...
	'Position',[205 088 083 024].*WS,...
	'ToolTipString','rendering utilities...',...
	'FontSize',FS(9),...
	'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{	'spm_render;',...
			'spm_surf;'	})

uicontrol(Fmenu,'Style','PopUp',...
	'String',Modalities,...
	'ToolTipString','change modality PET<->fMRI',...
	'Tag','Modality',...
	'Position',[298 088 082 024].*WS,...
	'CallBack',[...
		'if isempty(get(gco,''UserData'')) | ',...
		'get(gco,''Value'')~=get(gco,''UserData''),',...
			'spm(''ChMod'',get(gco,''Value'')),',...
		'end'],...
	'UserData','SPM_Modality',...
	'Interruptible','off')

%-Utility buttons (second line)
%-----------------------------------------------------------------------
%-Toolbox pulldown
xTB = spm('TBs');
if isempty(xTB)
    uicontrol(Fmenu,'String','Toolboxes...',...
	'Position',[020 054 082 024].*WS,...
	'ToolTipString','Additional Toolboxes',...
	'FontSize',FS(9),...
	'UserData','MenuWin_Toolboxes_PullDown',...
	'CallBack','','Enable','off')
else
    uicontrol(Fmenu,'Style','PopUp',...
	'String',{'Toolboxes...',xTB.name},...
       'ToolTipString','Additional Toolboxes',...
       'Position',[020 054 082 024].*WS,...
       'FontSize',FS(9),...
	'CallBack',...
        ['spm(''TBlaunch'',get(gcbo,''UserData''),get(gcbo,''Value'')-1), ',...
        	'set(gcbo,''Value'',1)'],...
       'UserData', xTB)
end

uicontrol(Fmenu,...
	'String','PPIs',...
	'Position',[112 054 083 024].*WS,...
	'ToolTipString',...
	'Forms PPI regressors for fMRI using hemodynamic deconvolution',...
	'FontSize',FS(9),...
	'UserData','spm_peb_ppi',...
	'CallBack','spm_peb_ppi')

uicontrol(Fmenu,'String','ImCalc',...
	'Position',[205 054 083 024].*WS,...
	'ToolTipString','image calculator',...
	'FontSize',FS(9),...
	'UserData','spm_imcalc_ui',...
	'CallBack','spm_imcalc_ui;')

uicontrol(Fmenu,'String','Bias cor',...
	'Position',[298 054 082 024].*WS,...
	'ToolTipString','Bias correction',...
	'FontSize',FS(9),...
	'UserData','spm_bias_ui',...
	'CallBack','spm_bias_ui;')

%-Utility buttons (third line)
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Help',...
	'Position',[020 020 082 024].*WS,...
	'ToolTipString','launch SPMhelp browser',...
	'CallBack','spm_help;',...
	'UserData','spm_help',...
	'ForeGroundColor','g')

uicontrol(Fmenu,'Style','PopUp',...
	'String','Utils...|CD|PWD|delete files|Show SPM|Run mFile|SPMweb',...
	'ToolTipString','misc SPM utilities',...
	'Position',[112 020 083 024].*WS,...
	'FontSize',FS(9),...
	'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{	[...
	    'spm(''FnBanner'',''CD'');',...
	       'cd(spm_get(-1,''*'',''Select new working directory'',pwd)),',...
	       'spm(''alert"'',',...
	    	    '{''New working directory:'',[''    '',pwd]},',...
	    	    '''CD'',sqrt(-1));'],...
	        ['spm(''alert"'',',...
	    	'{''Present working directory:'',[''    '',pwd]},',...
	    	'''PWD'',0);'],...
	    'spm(''GUI_FileDelete'')',...
	    'spm(''Show'');',...
	    'run(spm_get(1,''*.m'',''Select mFile to run''))',...
	    'web(''http://www.fil.ion.ucl.ac.uk/spm'')' } )

uicontrol(Fmenu,'String','Defaults',...
	'Position',[205 020 083 024].*WS,...
	'ToolTipString','adjust default SPM behavior for this session',...
	'FontSize',FS(9),...
	'CallBack','spm_defaults_edit;')

uicontrol(Fmenu,'String','Quit',...
	'Position',[298 020 082 024].*WS,...
	'ToolTipString','exit SPM',...
	'ForeGroundColor','r',...
	'Interruptible','off',...
	'UserData','QuitSPM',...
	'CallBack','spm(''Quit''), clear all')

%-----------------------------------------------------------------------
set(Fmenu,'CloseRequestFcn','spm(''Quit'')')
set(Fmenu,'Visible',Vis)
varargout = {Fmenu};



%=======================================================================
case 'createintwin'                      %-Create SPM interactive window
%=======================================================================
% Finter = spm('CreateIntWin',Vis)
if nargin<2, Vis='on'; else, Vis=varargin{2}; end

%-Close any existing 'Interactive' 'Tag'ged windows
delete(spm_figure('FindWin','Interactive'))

FS   = spm('FontSizes');			%-Scaled font sizes
PF   = spm_platform('fonts');			%-Font names (for this platform)
Rect = spm('WinSize','Interactive');		%-Interactive window rectangle

%-Create SPM Interactive window
Finter = figure('IntegerHandle','off',...
	'Tag','Interactive',...
	'Name','','NumberTitle','off',...
	'Position',Rect,...
	'Resize','off',...
	'Color',[1 1 1]*.7,...
	'MenuBar','none',...
	'DefaultTextFontName',PF.helvetica,...
	'DefaultTextFontSize',FS(10),...
	'DefaultAxesFontName',PF.helvetica,...
	'DefaultUicontrolBackgroundColor',[1 1 1]*.7,...
	'DefaultUicontrolFontName',PF.helvetica,...
	'DefaultUicontrolFontSize',FS(10),...
	'DefaultUicontrolInterruptible','on',...
	'Renderer', 'zbuffer',...
	'Visible',Vis);
varargout = {Finter};


%=======================================================================
case 'chmod'                            %-Change SPM modality PET<->fMRI
%=======================================================================
% spm('ChMod',Modality)

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality = ''; else, Modality = varargin{2}; end
[Modality,ModNum] = spm('CheckModality',Modality);

if strcmp(Modality,'PET'), OModality = 'FMRI'; else, OModality='PET'; end


%-Sort out global defaults
%-----------------------------------------------------------------------
spm('defaults',Modality)

%-Sort out visability of appropriate controls on Menu window
%-----------------------------------------------------------------------
Fmenu = spm_figure('FindWin','Menu');
if isempty(Fmenu), error('SPM Menu window not found'), end

set(findobj(Fmenu,'Tag',OModality),'Visible','off')
set(findobj(Fmenu,'Tag', Modality),'Visible','on')
set(findobj(Fmenu,'Tag','Modality'),'Value',ModNum,'UserData',ModNum)


%=======================================================================
case 'defaults'                 %-Set SPM defaults (as global variables)
%=======================================================================
% spm('defaults',Modality)
global defaults
if isempty(defaults), spm_defaults, end;

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else, Modality=varargin{2}; end
Modality          = spm('CheckModality',Modality);
defaults.modality = Modality;
defaults.SWD      = spm('Dir');				% SPM directory
defaults.TWD      = spm_platform('tempdir');		% Temp directory
	
%-Set Modality specific default (global variables)
%-----------------------------------------------------------------------
global UFp
if strcmp(lower(defaults.modality),'pet')
	UFp	= defaults.stats.pet.ufp;		% Upper tail F-prob
elseif strcmp(lower(defaults.modality),'fmri')
	UFp	= defaults.stats.fmri.ufp;		% Upper tail F-prob
elseif strcmp(lower(defaults.modality),'unknown')
else
	error('Illegal Modality')
end


%=======================================================================
case 'quit'                                      %-Quit SPM and clean up
%=======================================================================
% spm('Quit')
%-----------------------------------------------------------------------
delete(get(0,'Children'))
clc
fprintf('Bye...\n\n')


%=======================================================================
case 'checkmodality'              %-Check & canonicalise modality string
%=======================================================================
% [Modality,ModNum] = spm('CheckModality',Modality)
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else, Modality=upper(varargin{2}); end
if isempty(Modality)
	global defaults
	if isfield(defaults,'modality'); Modality = defaults.modality;
	else, Modality = 'UNKNOWN'; end
end
if isstr(Modality)
	ModNum = find(all(Modalities(:,1:length(Modality))'==...
			Modality'*ones(1,size(Modalities,1))));
else
	if ~any(Modality==[1:size(Modalities,1)])
		Modality = 'ERROR';
		ModNum   = [];
	else
		ModNum   = Modality;
		Modality = deblank(Modalities(ModNum,:));
	end
end

if isempty(ModNum), error('Unknown Modality'), end
varargout = {upper(Modality),ModNum};


%=======================================================================
case {'winscale','getwinscale'}  %-Window scale factors (to fit display)
%=======================================================================
% WS = spm('WinScale')
if strcmp(lower(Action),'getwinscale')
	warning('spm(''GetWinScale'' GrandFathered, use ''WinScale''')
end
S0   = get(0,'ScreenSize');
if all(S0==1), error('Can''t open any graphics windows...'), end
varargout = {[S0(3)/1152 (S0(4)-50)/900 S0(3)/1152 (S0(4)-50)/900]};


%=======================================================================
case {'fontsize','fontsizes','fontscale'}                 %-Font scaling
%=======================================================================
% [FS,sf] = spm('FontSize',FS)
% [FS,sf] = spm('FontSizes',FS)
% sf = spm('FontScale')
if nargin<3, c=0; else, c=1; end
if nargin<2, FS=[1:36]; else, FS=varargin{2}; end

sf  = 1 + 0.85*(min(spm('WinScale'))-1);

if strcmp(lower(Action),'fontscale')
	varargout = {sf};
else
	varargout = {ceil(FS*sf),sf};
end


%=======================================================================
case 'winsize'                 %-Standard SPM window locations and sizes
%=======================================================================
% Rect = spm('WinSize',Win,raw)
if nargin<3, raw=0; else, raw=1; end
if nargin<2, Win=''; else, Win=varargin{2}; end

Rect = [	[108 466 400 445];...
		[108 045 400 395];...
		[515 045 600 865] ];

WS = spm('WinScale');

if isempty(Win)
	WS = ones(3,1)*WS;
elseif upper(Win(1))=='M'
	%-Menu window
	Rect = Rect(1,:);
elseif upper(Win(1))=='I'
	%-Interactive window
	Rect = Rect(2,:);
elseif upper(Win(1))=='G'
	%-Graphics window
	Rect = Rect(3,:);
elseif Win(1)=='0'
	%-Root workspace
	Rect = get(0,'ScreenSize');
else
	error('Unknown Win type');
end

if ~raw, Rect = Rect.*WS; end
varargout = {Rect};


%=======================================================================
case 'dir'                           %-Identify specific (SPM) directory
%=======================================================================
% spm('Dir',Mfile)
%-----------------------------------------------------------------------
if nargin<2, Mfile='spm'; else, Mfile=varargin{2}; end
SPMdir = which(Mfile);
if isempty(SPMdir)			%-Not found or full pathname given
	if exist(Mfile,'file')==2	%-Full pathname
		SPMdir = Mfile;
	else
		error(['Can''t find ',Mfile,' on MATLABPATH']);
	end
end
varargout = {spm_str_manip(SPMdir,'H')};


%=======================================================================
case 'ver'                                                 %-SPM version
%=======================================================================
% SPMver = spm('Ver',Mfile,ReDo,Cache,Con)
if nargin<5, Con=[]; else, Con=varargin{5}; end
if isempty(Con), Con=1; end
if nargin<4, Cache=[]; else, Cache=varargin{4}; end
if isempty(Cache), Cache=1; end
if nargin<3, ReDo=[]; else, ReDo=varargin{3}; end
if isempty(ReDo), ReDo=0; end
if nargin<2, Mfile=''; else, Mfile=varargin{2}; end
if isempty(Mfile), Mfile='spm'; end

xVname = [upper(spm_str_manip(Mfile,'rt')),'_VER'];

%-See if version info exists in global variable
%-----------------------------------------------------------------------
xV = spm('GetGlobal',xVname);
if ~ReDo & ~isempty(xV)
	if isstruct(xV) & isfield(xV,'v') & isfield(xV,'c')
		varargout = {xV.v,xV.c};
		return
	end
end

%-Work version out from file
%-----------------------------------------------------------------------
if Con
	Vfile = fullfile(spm('Dir',Mfile),'Contents.m');
	skip = 0;	%-Don't skip first line
else
	Vfile = which(Mfile);
	if isempty(Vfile), error(['Can''t find ',Mfile,' on MATLABPATH']); end
	skip = 1;	%-Skip first line
end
if exist(Vfile)
	fid = fopen(Vfile,'r');
	str = fgets(fid);
	if skip, str=fgets(fid); end
	fclose(fid);
	str(1:max(1,min(find(str~='%' & str~=' '))-1))=[];
	tmp = min(find(str==10|str==32));
	v = str(1:tmp-1);
	if str(tmp)==32
		c = str(tmp+1:tmp+min(find(str(tmp+1:end)==10))-1);
	else
		c = '(c) Copyright reserved';
	end
else
	v = 'SPM';
	c = '(c) Copyright reserved';
end

%-Store version info in global variable
%-----------------------------------------------------------------------
if Cache
	eval(['global ',xVname])
	eval([xVname,' = struct(''v'',v,''c'',c);'])
end

varargout = {v,c};


%=======================================================================
case 'tbs'                                %-Identify installed toolboxes
%=======================================================================
% xTB = spm('TBs')

Tdir   = fullfile(spm('Dir'),'toolbox');

if exist(Tdir,'dir')
	[null,tmp] = spm_list_files(Tdir,'-');
	tmp = cellstr(tmp(tmp(:,1)~='.',:))';
else
	tmp = {};
end

tboxs = {};
tprog = {};
tdirs = {};
for tbox = tmp
	tbox = char(tbox);
	tdir = fullfile(Tdir,tbox);
	[fn,null]=spm_list_files(tdir,['*',tbox,'.m']);
	if size(fn,1)==1,
		tboxs = [tboxs, {strrep(tbox,'_','')}];
		tprog = [tprog, {spm_str_manip(fn,'r')}];
		tdirs = [tdirs, {tdir}];
	end
end

varargout = {struct('name',tboxs,'prog',tprog,'dir',tdirs)};


%=======================================================================
case 'tblaunch'                                  %-Launch an SPM toolbox
%=======================================================================
% xTB = spm('TBlaunch',xTB,i)

if nargin<3, i=1; else, i=varargin{3}; end
if nargin<2, error('insufficient arguments'), end

%-Addpath (& report)
%-----------------------------------------------------------------------
if i>0,
	if isempty(findstr(varargin{2}(i).dir,path))
		addpath(varargin{2}(i).dir,'-begin')
		spm('alert"',{'Toolbox directory prepended to MatLab path:',...
			varargin{2}(i).dir},...
			[varargin{2}(i).name,' toolbox'],1);
	end

	%-Launch
	%-----------------------------------------------------------------------
	evalin('base',varargin{2}(i).prog)
end;


%=======================================================================
case 'colour'                                     %-SPM interface colour
%=======================================================================
% spm('Colour')
%-----------------------------------------------------------------------
%-Pre-developmental livery
% varargout = {[1.0,0.2,0.3],'fightening red'};
%-Developmental livery
% varargout = {[0.7,1.0,0.7],'flourescent green'};
%-Alpha release livery
% varargout = {[0.9,0.9,0.5],'over-ripe banana'};
%-Beta release livery
  varargout = {[0.9 0.8 0.9],'blackcurrant purple'};
%-Distribution livery
% varargout = {[0.8 0.8 1.0],'vile violet'};;


%=======================================================================
case 'getglobal'                           %-Get global variable cleanly
%=======================================================================
% varargout = spm('GetGlobal',varargin)
wg = who('global');
for i=1:nargin-1
	if any(strcmp(wg,varargin{i+1}))
		eval(['global ',varargin{i+1},', tmp=',varargin{i+1},';'])
		varargout{i} = tmp;
	else
		varargout{i} = [];
	end
end

%=======================================================================
case {'cmdline','isgcmdline'}                   %-SPM command line mode?
%=======================================================================
% CmdLine = spm('CmdLine',CmdLine)
% isGCmdLine usage is Grandfathered
if nargin<2, CmdLine=[]; else, CmdLine = varargin{2}; end
if isempty(CmdLine),
	global defaults
	if ~isempty(defaults) & isfield(defaults,'cmdline'),
		CmdLine = defaults.cmdline;
	else,
		CmdLine = 0;
	end;
end
varargout = {CmdLine * (get(0,'ScreenDepth')>0)};

%=======================================================================
case 'mlver'                       %-MatLab major & point version number
%=======================================================================
% v = spm('MLver')
v = version; tmp = find(v=='.');
if length(tmp)>1, varargout={v(1:tmp(2)-1)}; end

%=======================================================================
case 'setcmdwinlabel'      %-Set command window label (Sun OpenWin only)
%=======================================================================
% spm('SetCmdWinLabel',WinStripe,IconLabel)
%-----------------------------------------------------------------------

%-Only label Sun command tools
%-----------------------------------------------------------------------
Term        = getenv('TERM');
if ~strcmp(Term,'sun-cmd'), return, end

%-Work out label text
%-----------------------------------------------------------------------
User        = spm('GetUser');
[null,Host] = unix('echo `hostname` | sed -e ''s/\..*$//''');
Host        = Host(1:length(Host)-1); 
v           = spm('MLver');

if nargin<3, IconLabel = ['MatLab',v(1)]; end
if nargin<2, WinStripe = [User,' - ',Host,' : MatLab ',v]; end

%-Set window stripe
%-----------------------------------------------------------------------
disp([']l' WinStripe '\]L' IconLabel '\'])


%=======================================================================
case 'popupcb'               %-Callback handling utility for PopUp menus
%=======================================================================
% spm('PopUpCB',h)
if nargin<2, h=gcbo; else, h=varargin{2}; end
v   = get(h,'Value');
if v==1, return, end
set(h,'Value',1)
CBs = get(h,'UserData');
evalin('base',CBs{v-1})


%=======================================================================
case 'getuser'                                           %-Get user name
%=======================================================================
% str = spm('GetUser',fmt)
str = spm_platform('user');
if ~isempty(str) & nargin>1, str = sprintf(varargin{2},str); end
varargout = {str};


%=======================================================================
case 'beep'                                %-Emit a keyboard "bell" beep
%=======================================================================
% spm('Beep')
fprintf('%c',7)


%=======================================================================
case 'time'                          %-Return formatted date/time string
%=======================================================================
% [timestr, date_vec] = spm('Time')
tmp = clock;
varargout = {sprintf('%02d:%02d:%02d - %02d/%02d/%4d',...
			tmp(4),tmp(5),floor(tmp(6)),tmp(3),tmp(2),tmp(1)),...
		tmp};


%=======================================================================
case 'pointer'                 %-Set mouse pointer in all MatLab windows
%=======================================================================
% spm('Pointer',Pointer)
if nargin<2, Pointer='Arrow'; else, Pointer=varargin{2}; end
set(get(0,'Children'),'Pointer',Pointer)


%=======================================================================
case {'alert','alert"','alert*','alert!'}                %-Alert dialogs
%=======================================================================
% h = spm('alert',Message,Title,CmdLine,wait)

%- Globals 
%-----------------------------------------------------------------------
if nargin<5, wait    = 0;  else, wait    = varargin{5}; end
if nargin<4, CmdLine = []; else, CmdLine = varargin{4}; end
if nargin<3, Title   = ''; else, Title   = varargin{3}; end
if nargin<2, Message = ''; else, Message = varargin{2}; end
Message = cellstr(Message);

if isreal(CmdLine)
	CmdLine  = spm('CmdLine',CmdLine);
	CmdLine2 = 0;
else
	CmdLine  = spm('CmdLine');
	CmdLine2 = 1;
end
timestr = spm('Time');
SPMv    = spm('ver');

switch(lower(Action))
case 'alert',	icon = 'none';	str = '--- ';
case 'alert"',	icon = 'help';	str = '~ - ';
case 'alert*',	icon = 'error'; str = '* - ';
case 'alert!',	icon = 'warn';	str = '! - ';
end

if CmdLine | CmdLine2
	Message(strcmp(Message,'')) = {' '};
	tmp = sprintf('%s: %s',SPMv,Title);
	fprintf('\n    %s%s  %s\n\n',str,tmp,repmat('-',1,62-length(tmp)))
	fprintf('        %s\n',Message{:})
	fprintf('\n        %s  %s\n\n',repmat('-',1,62-length(timestr)),timestr)
	h = [];
end

if ~CmdLine
	tmp = max(size(char(Message),2),42) - length(SPMv) - length(timestr);
	str = sprintf('%s  %s  %s',SPMv,repmat(' ',1,tmp-4),timestr);
	h   = msgbox([{''};Message(:);{''};{''};{str}],...
		sprintf('%s%s: %s',SPMv,spm('GetUser',' (%s)'),Title),...
		icon,'non-modal');
	drawnow
	set(h,'windowstyle','modal');
end

if wait
	if isempty(h)
		input('        press ENTER to continue...');
	else
		uiwait(h)
		h = [];
	end
end

if nargout, varargout = {h}; end


%=======================================================================
case {'fnbanner','sfnbanner','ssfnbanner'}  %-Text banners for functions
%=======================================================================
% SPMid = spm('FnBanner', Fn,FnV)
% SPMid = spm('SFnBanner',Fn,FnV)
% SPMid = spm('SSFnBanner',Fn,FnV)
time = spm('time');
str  = spm('ver');
if nargin>=2, str = [str,': ',varargin{2}]; end
if nargin>=3, str = [str,' (v',varargin{3},')']; end

switch lower(Action)
case 'fnbanner'
	tab = '';
	wid = 72;
	lch = '=';
case 'sfnbanner'
	tab = sprintf('\t');
	wid = 72-8;
	lch = '-';
case 'ssfnbanner'
	tab = sprintf('\t\t');
	wid = 72-2*8;
	lch = '-';
end

fprintf('\n%s%s',tab,str)
fprintf('%c',' '*ones(1,wid-length([str,time])))
fprintf('%s\n%s',time,tab)
fprintf('%c',lch*ones(1,wid)),fprintf('\n')
varargout = {str};


%=======================================================================
case 'fnuisetup'                %-Robust UI setup for main SPM functions
%=======================================================================
% [Finter,Fgraph,CmdLine] = spm('FnUIsetup',Iname,bGX,CmdLine)
if nargin<4, CmdLine=spm('CmdLine'); else, CmdLine=varargin{4}; end
if nargin<3, bGX=1; else, bGX=varargin{3}; end
if nargin<2, Iname=''; else, Iname=varargin{2}; end
if CmdLine
	Finter = spm_figure('FindWin','Interactive');
	if ~isempty(Finter), spm_figure('Clear',Finter), end
	%if ~isempty(Iname), fprintf('%s:\n',Iname), end
else
	Finter = spm_figure('GetWin','Interactive');
	spm_figure('Clear',Finter)
	if ~isempty(Iname)
		str = sprintf('%s (%s): %s',spm('ver'),spm('GetUser'),Iname);
	else
		str = '';
	end
	set(Finter,'Name',str)
end

if bGX
	Fgraph = spm_figure('GetWin','Graphics');
	spm_figure('Clear',Fgraph)
else
	Fgraph = spm_figure('FindWin','Graphics');
end
varargout = {Finter,Fgraph,CmdLine};	


%=======================================================================
case 'figname'                                %-Robust SPM figure naming
%=======================================================================
% F = spm('FigName',Iname,F,CmdLine)
if nargin<4, CmdLine=spm('CmdLine'); else, CmdLine=varargin{4}; end
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, Iname=''; else, Iname=varargin{2}; end

%if ~isempty(Iname), fprintf('\t%s\n',Iname), end
if CmdLine, varargout={[]}; return, end
F = spm_figure('FindWin',F);
if ~isempty(F) & ~isempty(Iname)
	set(F,'Name',sprintf('%s (%s): %s',spm('ver'),spm('GetUser'),Iname))
end
varargout={F};


%=======================================================================
case 'gui_filedelete'                                %-GUI file deletion
%=======================================================================
% spm('GUI_FileDelete')
P = spm_get(Inf,'*',{'Select file(s) to delete'},pwd,0);
n = length(P);
if n==0
	spm('alert"','Nothing selected to delete!','file delete',0);
	return
elseif n<4
	str=[{' '};P];
elseif n<11
	str=[{' '};P;{' ';sprintf('(%d files)',n)}];
else
	str=[{' '};P(1:min(n,10));{'...';' ';sprintf('(%d files)',n)}];
end
if spm_input(str,-1,'bd','delete|cancel',[1,0],[],'confirm file delete')
	spm_unlink(P{:})
	spm('alert"',P,'file delete',1);
end


%=======================================================================
case 'show'                   %-Bring visible MatLab windows to the fore
%=======================================================================
% Fs = spm('Show')
cF = get(0,'CurrentFigure');
Fs = get(0,'Children');
Fs = findobj(Fs,'flat','Visible','on');
for F=Fs', figure(F), end
set(0,'CurrentFigure',cF)
spm('FnBanner','GUI show');
varargout={Fs};


%=======================================================================
case 'clear'                                             %-Clear SPM GUI
%=======================================================================
% spm('Clear',Finter, Fgraph)
if nargin<3, Fgraph='Graphics'; else, Fgraph=varargin{3}; end
if nargin<2, Finter='Interactive'; else, Finter=varargin{2}; end
spm_figure('Clear',Fgraph)
spm_figure('Clear',Finter)
spm('Pointer','Arrow')
spm_get('Initialise','reset');
spm_conman('Initialise','reset');
clc, spm('FnBanner','GUI cleared');
fprintf('\n');
%evalin('Base','clear')


%=======================================================================
case 'help'                                  %-Pass through for spm_help
%=======================================================================
% spm('Help',varargin)
if nargin>1, spm_help(varargin{2:end}), else, spm_help, end


%=======================================================================
otherwise                                        %-Unknown action string
%=======================================================================
error('Unknown action string')

%=======================================================================
end
