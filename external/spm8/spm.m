function varargout=spm(varargin)
% SPM: Statistical Parametric Mapping (startup function)
%_______________________________________________________________________
%  ___  ____  __  __
% / __)(  _ \(  \/  )  
% \__ \ )___/ )    (   Statistical Parametric Mapping
% (___/(__)  (_/\/\_)  SPM - http://www.fil.ion.ucl.ac.uk/spm/
%_______________________________________________________________________
%
% SPM (Statistical Parametric Mapping) is a package for the analysis
% functional brain mapping experiments. It is the in-house package of
% the Wellcome Trust Centre for Neuroimaging, and is available to the
% scientific community as copyright freeware under the terms of the
% GNU General Public Licence.
% 
% Theoretical, computational and other details of the package are
% available in SPM's "Help" facility. This can be launched from the
% main SPM Menu window using the "Help" button, or directly from the
% command line using the command `spm_help`.
%
% Details of this release are available via the "About SPM" help topic
% (file spm.man), accessible from the SPM splash screen.  (Or type
% `spm_help spm.man` in the MATLAB command window)
% 
% This spm function initialises the default parameters, and displays a
% splash screen with buttons leading to the PET, fMRI and M/EEG
% modalities. Alternatively, `spm('pet')`, `spm('fmri')`, `spm('eeg')`
% (equivalently `spm pet`, `spm fmri` and `spm eeg`) lead directly to 
% the respective modality interfaces.
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
% (managed by spm_select). See the help on spm_input.m and spm_select.m for
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
% Wellcome Trust Centre for Neuroimaging

%-SVN ID and authorship of this program...
%-----------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id$


%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take. Recall   )
%( MATLAB's command-function duality: `spm Welcome` is equivalent to   )
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
% Prints ASCII welcome banner in MATLAB command window.
%
% FORMAT spm('PET') spm('FMRI') spm('EEG')
% Closes all windows and draws new Menu, Interactive, and Graphics
% windows for an SPM session. The buttons in the Menu window launch the
% main analysis routines.
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
% FORMAT [c,cName] = spm('Colour')
% Returns the RGB triple and a description for the current en-vogue SPM
% colour, the background colour for the Menu and Help windows.
%
% FORMAT F = spm('FigName',Iname,F,CmdLine)
% Set name of figure F to "SPMver (User): Iname" if ~CmdLine
% Robust to absence of figure.
% Iname      - Name for figure
% F (input)  - Handle (or 'Tag') of figure to name [default 'Interactive']
% CmdLine    - CommandLine usage? [default spm('CmdLine')]
% F (output) - Handle of figure named
%
% FORMAT Fs = spm('Show')
% Opens all SPM figure windows (with HandleVisibility) using `figure`.
%   Maintains current figure.
% Fs - vector containing all HandleVisible figures (i.e. get(0,'Children'))
%
% FORMAT spm('Clear',Finter, Fgraph)
% Clears and resets SPM-GUI, clears and timestamps MATLAB command window.
% Finter  - handle or 'Tag' of 'Interactive' figure [default 'Interactive']
% Fgraph  - handle or 'Tag' of 'Graphics' figure [default 'Graphics']
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
% FORMAT SPMdir = spm('Dir',Mfile)
% Returns the directory containing the version of spm in use,
% identified as the first in MATLABPATH containing the Mfile spm (this
% file) (or Mfile if specified).
%
% FORMAT [v,r] = spm('Ver',Mfile,ReDo)
% Returns the current version (v) and release (r) of file Mfile. This 
% corresponds to the Last changed Revision number extracted from the 
% Subversion Id tag.
% If Mfile is absent or empty then it returns the current SPM version (v)
% and release (r), extracted from the file Contents.m in the SPM directory
% (these information are cached in a persistent variable to enable repeat
% use without recomputation).
% If Redo [default false] is true, then the cached current SPM information
% are not used but recomputed (and recached).
%
% FORMAT v = spm('MLver')
% Returns MATLAB version, truncated to major & minor revision numbers
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
% FORMAT spm('PopUpCB',h)
% Callback handler for PopUp UI menus with multiple callbacks as cellstr UserData
%
% FORMAT str = spm('GetUser',fmt)
% Returns current users login name, extracted from the hosting environment
% fmt   - format string: If USER is defined then sprintf(fmt,USER) is returned
%
% FORMAT spm('Beep')
% Plays the keyboard beep!
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
% FORMAT spm('GUI_FileDelete')
% CallBack for GUI for file deletion, using spm_select and confirmation dialogs
%
% FORMAT spm('Clean')
% Clear all variables, globals, functions, MEX links and class definitions.
%
% FORMAT spm('Help',varargin)
% Merely a gateway to spm_help(varargin) - so you can type "spm help"
%
% FORMAT spm('Quit')
% Quit SPM, delete all windows and clear the command window.
% 
%_______________________________________________________________________


%-Parameters
%-----------------------------------------------------------------------
Modalities = {'PET','FMRI','EEG'};

%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0, Action = 'Welcome'; else Action = varargin{1}; end


%=======================================================================
switch lower(Action), case 'welcome'             %-Welcome splash screen
%=======================================================================
% spm('Welcome')
spm_check_installation('basic');
defaults = spm('GetGlobal','defaults');
if isfield(defaults,'modality')
    spm(defaults.modality);
    return
end

%-Open startup window, set window defaults
%-----------------------------------------------------------------------
Fwelcome = openfig(fullfile(spm('Dir'),'spm_Welcome.fig'),'new','invisible');
set(Fwelcome,'name',sprintf('%s%s',spm('ver'),spm('GetUser',' (%s)')));
set(get(findobj(Fwelcome,'Type','axes'),'children'),'FontName',spm_platform('Font','Times'));
set(findobj(Fwelcome,'Tag','SPM_VER'),'String',spm('Ver'));
RectW = spm('WinSize','W',1); Rect0 = spm('WinSize','0',1);
set(Fwelcome,'Units','pixels', 'Position',...
    [Rect0(1)+(Rect0(3)-RectW(3))/2, Rect0(2)+(Rect0(4)-RectW(4))/2, RectW(3), RectW(4)]);
set(Fwelcome,'Visible','on');

%=======================================================================
case 'asciiwelcome'                           %-ASCII SPM banner welcome
%=======================================================================
% spm('AsciiWelcome')
disp( ' ___  ____  __  __                                            ');
disp( '/ __)(  _ \(  \/  )                                           ');
disp( '\__ \ )___/ )    (   Statistical Parametric Mapping           ');
disp(['(___/(__)  (_/\/\_)  ',spm('Ver'),' - http://www.fil.ion.ucl.ac.uk/spm/']);
fprintf('\n');


%=======================================================================
case lower(Modalities)       %-Initialise SPM in PET, fMRI, EEG modality
%=======================================================================
% spm(Modality)
spm_check_installation('basic');
try, feature('JavaFigures',0); end

%-Initialisation and workspace canonicalisation
%-----------------------------------------------------------------------
local_clc;
spm('AsciiWelcome');                    fprintf('\n\nInitialising SPM');
Modality = upper(Action);                                  fprintf('.');
delete(get(0,'Children'));                                 fprintf('.');

%-Load startup global defaults
%-----------------------------------------------------------------------
spm_defaults;                                              fprintf('.');

%-Setup for batch system
%-----------------------------------------------------------------------
spm_jobman('initcfg');
spm_select('prevdirs',[spm('Dir') filesep]);

%-Draw SPM windows
%-----------------------------------------------------------------------
if ~spm('CmdLine')
    Fmenu  = spm('CreateMenuWin','off');                   fprintf('.');
    Finter = spm('CreateIntWin','off');                    fprintf('.');
else
    Fmenu  = [];
    Finter = [];
end
Fgraph = spm_figure('Create','Graphics','Graphics','off'); fprintf('.');
   
spm_figure('WaterMark',Finter,spm('Ver'),'',45);           fprintf('.');

Fmotd  = fullfile(spm('Dir'),'spm_motd.man');
if exist(Fmotd,'file'), spm_help('!Disp',Fmotd,'',Fgraph,spm('Ver')); end
                                                           fprintf('.');

%-Setup for current modality
%-----------------------------------------------------------------------
spm('ChMod',Modality);                                     fprintf('.');

%-Reveal windows
%-----------------------------------------------------------------------
set([Fmenu,Finter,Fgraph],'Visible','on');          fprintf('done\n\n');

%-Print present working directory
%-----------------------------------------------------------------------
fprintf('SPM present working directory:\n\t%s\n',pwd)


%=======================================================================
case 'chmod'                      %-Change SPM modality PET<->fMRI<->EEG
%=======================================================================
% spm('ChMod',Modality)
%-----------------------------------------------------------------------

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<2, Modality = ''; else Modality = varargin{2}; end
[Modality,ModNum] = spm('CheckModality',Modality);

%-Sort out global defaults
%-----------------------------------------------------------------------
spm('defaults',Modality);

%-Sort out visability of appropriate controls on Menu window
%-----------------------------------------------------------------------
Fmenu = spm_figure('FindWin','Menu');
if ~isempty(Fmenu)

    if strcmpi(Modality,'PET')
        set(findobj(Fmenu, 'Tag', 'FMRI'),    'Visible', 'off');
        set(findobj(Fmenu, 'Tag', 'EEG'),     'Visible', 'off');
        set(findobj(Fmenu, 'Tag', 'PETFMRI'), 'Visible', 'on' );
        set(findobj(Fmenu, 'Tag', 'PET'),     'Visible', 'on' );
    elseif strcmpi(Modality,'FMRI')
        set(findobj(Fmenu, 'Tag', 'EEG'),     'Visible', 'off');
        set(findobj(Fmenu, 'Tag', 'PET'),     'Visible', 'off');
        set(findobj(Fmenu, 'Tag', 'PETFMRI'), 'Visible', 'on' );
        set(findobj(Fmenu, 'Tag', 'FMRI'),    'Visible', 'on' );
    else
        set(findobj(Fmenu, 'Tag', 'PETFMRI'), 'Visible', 'off');
        set(findobj(Fmenu, 'Tag', 'PET'),     'Visible', 'off');
        set(findobj(Fmenu, 'Tag', 'FMRI'),    'Visible', 'off');
        set(findobj(Fmenu, 'Tag', 'EEG'),     'Visible', 'on' );
    end
    set(findobj(Fmenu,'Tag','Modality'),'Value',ModNum,'UserData',ModNum);
else
    warning('SPM Menu window not found');
end

%=======================================================================
case 'defaults'                  %-Set SPM defaults (as global variable)
%=======================================================================
% spm('defaults',Modality)
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else Modality=varargin{2}; end
Modality = spm('CheckModality',Modality);

%-Re-initialise, load defaults (from spm_defaults.m) and store modality
%-----------------------------------------------------------------------
clear global defaults
spm_get_defaults('modality',Modality);

%-Addpath modality-specific toolboxes
%-----------------------------------------------------------------------
if strcmpi(Modality,'EEG') && ~isdeployed
    addpath(fullfile(spm('Dir'),'external','fieldtrip'));
    ft_defaults;
    addpath(fullfile(spm('Dir'),'external','bemcp'));
    addpath(fullfile(spm('Dir'),'external','ctf'));
    addpath(fullfile(spm('Dir'),'external','eeprobe'));
    addpath(fullfile(spm('Dir'),'toolbox', 'dcm_meeg'));
    addpath(fullfile(spm('Dir'),'toolbox', 'spectral'));
    addpath(fullfile(spm('Dir'),'toolbox', 'Neural_Models'));
    addpath(fullfile(spm('Dir'),'toolbox', 'Beamforming'));
    addpath(fullfile(spm('Dir'),'toolbox', 'MEEGtools'));
end

%-Return defaults variable if asked
%-----------------------------------------------------------------------
if nargout, varargout = {spm_get_defaults}; end


%=======================================================================
case 'checkmodality'              %-Check & canonicalise modality string
%=======================================================================
% [Modality,ModNum] = spm('CheckModality',Modality)
%-----------------------------------------------------------------------
if nargin<2, Modality=''; else Modality=upper(varargin{2}); end
if isempty(Modality)
    try
        Modality = spm_get_defaults('modality');
    end
end
if ischar(Modality)
    ModNum = find(ismember(Modalities,Modality));
else
    if ~any(Modality == 1:length(Modalities))
        Modality = 'ERROR';
        ModNum   = [];
    else
        ModNum   = Modality;
        Modality = Modalities{ModNum};
    end
end
if isempty(ModNum)
    if isempty(Modality)
        fprintf('Modality is not set: use spm(''defaults'',''MOD''); ');
        fprintf('where MOD is one of PET, FMRI, EEG.\n');
    end
    error('Unknown Modality.');
end
varargout = {upper(Modality),ModNum};


%=======================================================================
case 'createmenuwin'                            %-Create SPM menu window
%=======================================================================
% Fmenu = spm('CreateMenuWin',Vis)
%-----------------------------------------------------------------------
if nargin<2, Vis='on'; else Vis=varargin{2}; end

%-Close any existing 'Menu' 'Tag'ged windows
%-----------------------------------------------------------------------
delete(spm_figure('FindWin','Menu'))
Fmenu = openfig(fullfile(spm('Dir'),'spm_Menu.fig'),'new','invisible');
set(Fmenu,'name',sprintf('%s%s: Menu',spm('ver'),spm('GetUser',' (%s)')));
S0 = spm('WinSize','0',1);
SM = spm('WinSize','M');
set(Fmenu,'Units','pixels', 'Position',[S0(1) S0(2) 0 0] + SM);

%-Set SPM colour
%-----------------------------------------------------------------------
set(findobj(Fmenu,'Tag', 'frame'),'backgroundColor',spm('colour'));
try
    if ismac
        b = findobj(Fmenu,'Style','pushbutton');
        set(b,'backgroundColor',get(b(1),'backgroundColor')+0.002);
    end
end

%-Set toolbox
%-----------------------------------------------------------------------
xTB       = spm('tbs');
if ~isempty(xTB)
    set(findobj(Fmenu,'Tag', 'Toolbox'),'String',{'Toolbox:' xTB.name });
    set(findobj(Fmenu,'Tag', 'Toolbox'),'UserData',xTB);
else
    set(findobj(Fmenu,'Tag', 'Toolbox'),'Visible','off')
end
set(Fmenu,'Visible',Vis);
varargout = {Fmenu};


%=======================================================================
case 'createintwin'                      %-Create SPM interactive window
%=======================================================================
% Finter = spm('CreateIntWin',Vis)
%-----------------------------------------------------------------------
if nargin<2, Vis='on'; else Vis=varargin{2}; end

%-Close any existing 'Interactive' 'Tag'ged windows
%-----------------------------------------------------------------------
delete(spm_figure('FindWin','Interactive'))
Finter = openfig(fullfile(spm('Dir'),'spm_Interactive.fig'),'new','invisible');
set(Finter,'name',spm('Ver'));
S0 = spm('WinSize','0',1);
SI = spm('WinSize','I');
set(Finter,'Units','pixels', 'Position',[S0(1) S0(2) 0 0] +  SI);
set(Finter,'Visible',Vis);
varargout = {Finter};


%=======================================================================
case 'fnuisetup'                %-Robust UI setup for main SPM functions
%=======================================================================
% [Finter,Fgraph,CmdLine] = spm('FnUIsetup',Iname,bGX,CmdLine)
%-----------------------------------------------------------------------
if nargin<4, CmdLine=spm('CmdLine'); else CmdLine=varargin{4}; end
if nargin<3, bGX=1; else bGX=varargin{3}; end
if nargin<2, Iname=''; else Iname=varargin{2}; end
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
case 'winscale'                  %-Window scale factors (to fit display)
%=======================================================================
% WS = spm('WinScale')
%-----------------------------------------------------------------------
S0 = spm('WinSize','0',1);

if all(ismember(S0(:),[0 1]))
    varargout = {[1 1 1 1]};
    return;
end

tmp   = [S0(3)/1152 (S0(4)-50)/900];

varargout = {min(tmp)*[1 1 1 1]};

% Make sure that aspect ratio is about right - for funny shaped screens
% varargout = {[S0(3)/1152 (S0(4)-50)/900 S0(3)/1152 (S0(4)-50)/900]};


%=======================================================================
case {'fontsize','fontsizes','fontscale'}                 %-Font scaling
%=======================================================================
% [FS,sf] = spm('FontSize',FS)
% [FS,sf] = spm('FontSizes',FS)
% sf = spm('FontScale')
%-----------------------------------------------------------------------
if nargin<2, FS=1:36; else FS=varargin{2}; end

offset     = 1;
%try, if ismac, offset = 1.4; end; end

sf  = offset + 0.85*(min(spm('WinScale'))-1);

if strcmpi(Action,'fontscale')
    varargout = {sf};
else
    varargout = {ceil(FS*sf),sf};
end


%=======================================================================
case 'winsize'                 %-Standard SPM window locations and sizes
%=======================================================================
% Rect = spm('WinSize',Win,raw)
%-----------------------------------------------------------------------
if nargin<3, raw=0; else raw=1; end
if nargin<2, Win=''; else Win=varargin{2}; end

Rect = [[108 466 400 445];...
        [108 045 400 395];...
        [515 015 600 865];...
        [326 310 500 280]];

if isempty(Win)
    %-All windows
elseif upper(Win(1))=='M'
    %-Menu window
    Rect = Rect(1,:);
elseif upper(Win(1))=='I'
    %-Interactive window
    Rect = Rect(2,:);
elseif upper(Win(1))=='G'
    %-Graphics window
    Rect = Rect(3,:);
elseif upper(Win(1))=='W'
    %-Welcome window
    Rect = Rect(4,:);
elseif Win(1)=='0'
    %-Root workspace
    Rect = get(0, 'MonitorPosition');
    if all(ismember(Rect(:),[0 1]))
        warning('SPM:noDisplay','Unable to open display.');
    end
    if size(Rect,1) > 1 % Multiple Monitors
        %-Use Monitor containing the Pointer
        pl = get(0,'PointerLocation');
        Rect(:,[3 4]) = Rect(:,[3 4]) + Rect(:,[1 2]);
        w  = find(pl(1)>=Rect(:,1) & pl(1)<=Rect(:,3) &...
                  pl(2)>=Rect(:,2) & pl(2)<=Rect(:,4));
        if numel(w)~=1, w = 1; end
        Rect = Rect(w,:);
        %-Make sure that the format is [x y width height]
        Rect(1,[3 4]) = Rect(1,[3 4]) - Rect(1,[1 2]) + 1;
    end
else
    error('Unknown Win type');
end

if ~raw
    WS = repmat(spm('WinScale'),size(Rect,1),1);
    Rect = Rect.*WS; 
end

varargout = {Rect};


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
% varargout = {[0.9 0.8 0.9],'blackcurrant purple'};
%-Distribution livery
  varargout = {[0.8 0.8 1.0],'vile violet'};

try
    varargout = {spm_get_defaults('ui.colour'),'bluish'};
end
    

%=======================================================================
case 'figname'                                %-Robust SPM figure naming
%=======================================================================
% F = spm('FigName',Iname,F,CmdLine)
%-----------------------------------------------------------------------
if nargin<4, CmdLine=spm('CmdLine'); else CmdLine=varargin{4}; end
if nargin<3, F='Interactive'; else F=varargin{3}; end
if nargin<2, Iname=''; else Iname=varargin{2}; end

%if ~isempty(Iname), fprintf('\t%s\n',Iname), end
if CmdLine, varargout={[]}; return, end
F = spm_figure('FindWin',F);
if ~isempty(F) && ~isempty(Iname)
    set(F,'Name',sprintf('%s (%s): %s',spm('ver'),spm('GetUser'),Iname))
end
varargout={F};


%=======================================================================
case 'show'                   %-Bring visible MATLAB windows to the fore
%=======================================================================
% Fs = spm('Show')
%-----------------------------------------------------------------------
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
%-----------------------------------------------------------------------
if nargin<3, Fgraph='Graphics'; else Fgraph=varargin{3}; end
if nargin<2, Finter='Interactive'; else Finter=varargin{2}; end
spm_figure('Clear',Fgraph)
spm_figure('Clear',Finter)
spm('Pointer','Arrow')
spm_conman('Initialise','reset');
local_clc;
fprintf('\n');
%evalin('base','clear')


%=======================================================================
case {'fnbanner','sfnbanner','ssfnbanner'}  %-Text banners for functions
%=======================================================================
% SPMid = spm('FnBanner', Fn,FnV)
% SPMid = spm('SFnBanner',Fn,FnV)
% SPMid = spm('SSFnBanner',Fn,FnV)
%-----------------------------------------------------------------------
time = spm('time');
str  = spm('ver');
if nargin>=2, str = [str,': ',varargin{2}]; end
if nargin>=3 
    v = regexp(varargin{3},'\$Rev: (\d*) \$','tokens','once');
    if ~isempty(v)
        str = [str,' (v',v{1},')'];
    else
        str = [str,' (v',varargin{3},')'];
    end
end

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
fprintf('%c',repmat(' ',1,wid-length([str,time])))
fprintf('%s\n%s',time,tab)
fprintf('%c',repmat(lch,1,wid)),fprintf('\n')
varargout = {str};


%=======================================================================
case 'dir'                           %-Identify specific (SPM) directory
%=======================================================================
% spm('Dir',Mfile)
%-----------------------------------------------------------------------
if nargin<2, Mfile='spm'; else Mfile=varargin{2}; end
SPMdir = which(Mfile);
if isempty(SPMdir)             %-Not found or full pathname given
    if exist(Mfile,'file')==2  %-Full pathname
        SPMdir = Mfile;
    else
        error(['Can''t find ',Mfile,' on MATLABPATH']);
    end
end
SPMdir = fileparts(SPMdir);

% if isdeployed
%     ind = findstr(SPMdir,'_mcr')-1;
%     if ~isempty(ind)
%         % MATLAB 2008a/2009a doesn't need this
%         SPMdir = fileparts(SPMdir(1:ind(1)));
%     end
% end
varargout = {SPMdir};


%=======================================================================
case 'ver'                                                 %-SPM version
%=======================================================================
% [SPMver, SPMrel] = spm('Ver',Mfile,ReDo)
%-----------------------------------------------------------------------
if nargin > 3, warning('This usage of "spm ver" is now deprecated.'); end
if nargin ~= 3, ReDo = false; else ReDo = logical(varargin{3}); end
if nargin == 1 || (nargin > 1 && isempty(varargin{2}))
    Mfile = ''; 
else
    Mfile = which(varargin{2});
    if isempty(Mfile)
        error('Can''t find %s on MATLABPATH.',varargin{2});
    end
end

v = spm_version(ReDo);

if isempty(Mfile)
    varargout = {v.Release v.Version};
else
    unknown = struct('file',Mfile,'id','???','date','','author','');
    if ~isdeployed
        fp  = fopen(Mfile,'rt');
        if fp == -1, error('Can''t read %s.',Mfile); end
        str = fread(fp,Inf,'*uchar');
        fclose(fp);
        str = char(str(:)');
        r = regexp(str,['\$Id: (?<file>\S+) (?<id>[0-9]+) (?<date>\S+) ' ...
            '(\S+Z) (?<author>\S+) \$'],'names','once');
        if isempty(r), r = unknown; end
    else
        r = unknown;
    end
    varargout = {r(1).id v.Release};
end


%=======================================================================
case 'mlver'                       %-MATLAB major & point version number
%=======================================================================
% v = spm('MLver')
%-----------------------------------------------------------------------
v = version; tmp = find(v=='.');
if length(tmp)>1, varargout={v(1:tmp(2)-1)}; end


%=======================================================================
case 'tbs'                                %-Identify installed toolboxes
%=======================================================================
% xTB = spm('TBs')
%-----------------------------------------------------------------------

% Toolbox directory
%-----------------------------------------------------------------------
Tdir  = fullfile(spm('Dir'),'toolbox');

%-List of potential installed toolboxes directories
%-----------------------------------------------------------------------
if exist(Tdir,'dir')
    d = dir(Tdir); 
    d = {d([d.isdir]).name};
    d = {d{cellfun('isempty',regexp(d,'^\.'))}};
else
    d = {};
end


%-Look for a "main" M-file in each potential directory
%-----------------------------------------------------------------------
xTB = [];
for i = 1:length(d)
    tdir = fullfile(Tdir,d{i});
    fn   = cellstr(spm_select('List',tdir,['^.*' d{i} '\.m$']));

    if ~isempty(fn{1}),
        xTB(end+1).name = strrep(d{i},'_','');
        xTB(end).prog   = spm_str_manip(fn{1},'r');
        xTB(end).dir    = tdir;
    end

end

varargout{1} = xTB;


%=======================================================================
case 'tblaunch'                                  %-Launch an SPM toolbox
%=======================================================================
% xTB = spm('TBlaunch',xTB,i)
%-----------------------------------------------------------------------
if nargin < 3, i   = 1;          else i   = varargin{3}; end
if nargin < 2, xTB = spm('TBs'); else xTB = varargin{2}; end

if i > 0
    %-Addpath (& report)
    %-------------------------------------------------------------------
    if isempty(findstr(xTB(i).dir,path))
        if ~isdeployed, addpath(xTB(i).dir,'-begin'); end
        spm('alert"',{'Toolbox directory prepended to Matlab path:',...
            xTB(i).dir},...
            [xTB(i).name,' toolbox'],1);
    end

    %-Launch
    %-------------------------------------------------------------------
    evalin('base',xTB(i).prog);
end


%=======================================================================
case 'getglobal'                           %-Get global variable cleanly
%=======================================================================
% varargout = spm('GetGlobal',varargin)
%-----------------------------------------------------------------------
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
case {'cmdline'}                                %-SPM command line mode?
%=======================================================================
% CmdLine = spm('CmdLine',CmdLine)
%-----------------------------------------------------------------------
if nargin<2, CmdLine=[]; else CmdLine=varargin{2}; end
if isempty(CmdLine)
    defaults = spm('getglobal','defaults');
    if isfield(defaults,'cmdline')
        CmdLine = defaults.cmdline;
    else
        CmdLine = 0;
    end
end
varargout = {CmdLine || (get(0,'ScreenDepth')==0)};


%=======================================================================
case 'popupcb'               %-Callback handling utility for PopUp menus
%=======================================================================
% spm('PopUpCB',h)
%-----------------------------------------------------------------------
if nargin<2, h=gcbo; else h=varargin{2}; end
v   = get(h,'Value');
if v==1, return, end
set(h,'Value',1)
CBs = get(h,'UserData');
evalin('base',CBs{v-1})


%=======================================================================
case 'getuser'                                           %-Get user name
%=======================================================================
% str = spm('GetUser',fmt)
%-----------------------------------------------------------------------
str = spm_platform('user');
if ~isempty(str) && nargin>1, str = sprintf(varargin{2},str); end
varargout = {str};


%=======================================================================
case 'beep'                                         %-Produce beep sound
%=======================================================================
% spm('Beep')
%-----------------------------------------------------------------------
beep;


%=======================================================================
case 'time'                          %-Return formatted date/time string
%=======================================================================
% [timestr, date_vec] = spm('Time')
%-----------------------------------------------------------------------
tmp = clock;
varargout = {sprintf('%02d:%02d:%02d - %02d/%02d/%4d',...
             tmp(4),tmp(5),floor(tmp(6)),tmp(3),tmp(2),tmp(1)), tmp};


%=======================================================================
case 'memory'
%=======================================================================
% m = spm('Memory')
%-----------------------------------------------------------------------
maxmemdef = 200*1024*1024;
m         = maxmemdef;
varargout = {m};

%=======================================================================
case 'pointer'                 %-Set mouse pointer in all MATLAB windows
%=======================================================================
% spm('Pointer',Pointer)
%-----------------------------------------------------------------------
if nargin<2, Pointer='Arrow'; else  Pointer=varargin{2}; end
set(get(0,'Children'),'Pointer',Pointer)




%=======================================================================
case {'alert','alert"','alert*','alert!'}                %-Alert dialogs
%=======================================================================
% h = spm('alert',Message,Title,CmdLine,wait)
%-----------------------------------------------------------------------

%- Globals 
%-----------------------------------------------------------------------
if nargin<5, wait    = 0;  else wait    = varargin{5}; end
if nargin<4, CmdLine = []; else CmdLine = varargin{4}; end
if nargin<3, Title   = ''; else Title   = varargin{3}; end
if nargin<2, Message = ''; else Message = varargin{2}; end
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
    case 'alert',  icon = 'none';  str = '--- ';
    case 'alert"', icon = 'help';  str = '~ - ';
    case 'alert*', icon = 'error'; str = '* - ';
    case 'alert!', icon = 'warn';  str = '! - ';
end

if CmdLine || CmdLine2
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
case 'gui_filedelete'                                %-GUI file deletion
%=======================================================================
% spm('GUI_FileDelete')
%-----------------------------------------------------------------------
[P, sts] = spm_select(Inf,'.*','Select file(s) to delete');
if ~sts, return; end
P = cellstr(P);
n = numel(P);
if n==0 || (n==1 && isempty(P{1}))
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
    feval(@spm_unlink,P{:});
    spm('alert"',P,'file delete',1);
end


%=======================================================================
case 'clean'                                    %-Clean MATLAB workspace
%=======================================================================
% spm('Clean')
%-----------------------------------------------------------------------
evalin('base','clear all');
evalc('clear classes');


%=======================================================================
case 'help'                                  %-Pass through for spm_help
%=======================================================================
% spm('Help',varargin)
%-----------------------------------------------------------------------
if nargin>1, spm_help(varargin{2:end}), else spm_help, end


%=======================================================================
case 'quit'                                      %-Quit SPM and clean up
%=======================================================================
% spm('Quit')
%-----------------------------------------------------------------------
delete(get(0,'Children'));
local_clc;
fprintf('Bye for now...\n\n');


%=======================================================================
otherwise                                        %-Unknown action string
%=======================================================================
error('Unknown action string');

%=======================================================================
end


%=======================================================================
function local_clc                                %-Clear command window
%=======================================================================
if ~isdeployed
    clc;
end


%=======================================================================
function v = spm_version(ReDo)                    %-Retrieve SPM version
%=======================================================================
persistent SPM_VER;
v = SPM_VER;
if isempty(SPM_VER) || (nargin > 0 && ReDo)
    v = struct('Name','','Version','','Release','','Date','');
    try
        fid = fopen(fullfile(spm('Dir'),'Contents.m'),'rt');
        if fid == -1, error(str); end
        l1 = fgetl(fid); l2 = fgetl(fid);
        fclose(fid);
        l1 = strtrim(l1(2:end)); l2 = strtrim(l2(2:end));
        t  = textscan(l2,'%s','delimiter',' '); t = t{1};
        v.Name = l1; v.Date = t{4};
        v.Version = t{2}; v.Release = t{3}(2:end-1);
    catch
        if isdeployed
            % in deployed mode, M-files are encrypted
            % (but first two lines of Contents.m should be preserved)
            v.Name    = 'Statistical Parametric Mapping';
            v.Version = '8';
            v.Release = 'SPM8';
            v.Date    = date;
        else
            error('Can''t obtain SPM Revision information.');
        end
    end
    SPM_VER = v;
end
