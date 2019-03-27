function varargout=spm_figure(varargin)
% Setup and callback functions for Graphics window
% FORMAT varargout=spm_figure(varargin)
%
% spm_figure provides utility routines for using the SPM Graphics 
% interface. Most used syntaxes are listed here, see the embedded callback
% reference in the main body of this function, below the help text.
%
% FORMAT F = spm_figure('Create',Tag,Name,Visible)
% FORMAT F = spm_figure('FindWin',Tag)
% FORMAT F = spm_figure('GetWin',Tag)
% FORMAT spm_figure('Select',F)
% FORMAT spm_figure('Focus',F)
% FORMAT spm_figure('Clear',F,Tags)
% FORMAT spm_figure('Close',F)
% FORMAT spm_figure('Print',F)
% FORMAT spm_figure('WaterMark',F,str,Tag,Angle,Perm)
%
% FORMAT spm_figure('NewPage',hPage)
% FORMAT spm_figure('TurnPage',move,F)
% FORMAT spm_figure('DeletePageControls',F)
% FORMAT n = spm_figure('#page')
% FORMAT n = spm_figure('CurrentPage')
%__________________________________________________________________________
%
% spm_figure creates and manages the 'Graphics' window. This window and
% these facilities may be used independently of SPM, and any number of
% Graphics windows my be used within the same MATLAB session. (Though
% only one SPM 'Graphics' 'Tag'ed window is permitted).
%
% The Graphics window is provided with a menu bar at the top that
% facilitates editing and printing of the current graphic display.
% (This menu is also provided as a figure background "ContextMenu" - 
% right-clicking on the figure background should bring up the menu).
%
% "Print": Graphics windows with multi-page axes are printed page by page.
%
% "Clear": Clears the Graphics window. If in SPM usage (figure 'Tag'ed as
% 'Graphics') then all SPM windows are cleared and reset.
%
% "Colours":
% * gray, hot, pink, jet: Sets the colormap to selected item.
% * gray-hot, etc: Creates a 'split' colormap {128 x 3 matrix}.
%      The lower half is a gray scale and the upper half is selected
%      colormap  This colormap is used for viewing 'rendered' SPMs on a 
%      PET, MRI or other background images.
% Colormap effects:
% * Invert: Inverts (flips) the current color map.
% * Brighten and Darken: Brighten and Darken the current colourmap
%      using the MATLAB BRIGHTEN command, with  beta's of +0.2 and -0.2
%      respectively.
%
% For SPM usage, the figure should be 'Tag'ed as 'Graphics'.
%
% See also: spm_print, spm_clf
%__________________________________________________________________________
% Copyright (C) 1994-2015 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_figure.m 6409 2015-04-16 16:19:38Z guillaume $


%==========================================================================
% - FORMAT specifications for embedded CallBack functions
%==========================================================================
%
% FORMAT F = spm_figure
% [ShortCut] Defaults to Action 'Create'
%
% FORMAT F = spm_figure(F) - numeric F
% [ShortCut] Defaults to spm_figure('CreateBar',F)
%
% FORMAT F = spm_figure('Create',Tag,Name,Visible)
% Create a full length WhiteBg figure 'Tag'ed Tag (if specified),
% with a ToolBar and background context menu.
% Equivalent to spm_figure('CreateWin','Tag') and spm_figure('CreateBar')
% Tag     - 'Tag' string for figure.
% Name    - Name for window
% Visible - 'on' or 'off'
% F   - Figure used
%
% FORMAT F = spm_figure('FindWin',F)
% Finds window with 'Tag' or figure numnber F - returns empty F if not found
% F - (Input)  Figure to use [Optional] - 'Tag' string or figure number.
%   - Defaults to 'Graphics'
% F - (Output) Figure number (if found) or empty (if not).
%
% FORMAT F = spm_figure('GetWin',Tag)
% Like spm_figure('FindWin',Tag), except that if no such 'Tag'ged figure
% is found and 'Tag' is recognized, one is created. Further, the "got" 
% window is made current.
% Tag   - Figure 'Tag' to get, defaults to 'Graphics'
% F - Figure number (if found/created) or empty (if not).
%
% FORMAT spm_figure('Select',F)
% Set figure F as the current figure, so that any subsequent graphics
% commands are directed to it; however, this does not focus or raise the
% figure, so is appropriate for use on each iteration of a repetitive task.
%
% FORMAT spm_figure('Focus',F)
% Set figure F as the current figure and give it the focus (which depending
% on Operating System / Window Manager behavior might restore or raise the
% figure). This is equivalent to the usual figure(F) command, but in many
% cases spm_figure('Select',F) will be a more appropriate alternative.
%
% FORMAT spm_figure('Clear',F,Tags)
% Clears figure, leaving ToolBar (& other objects with invisible handles)
% Optional third argument specifies 'Tag's of objects to delete.
% If figure F is 'Tag'ged 'Interactive' (SPM usage), then the window
% name and pointer are reset.
% F - 'Tag' string or figure number of figure to clear, defaults to gcf
% Tags  - 'Tag's (string matrix or cell array of strings) of objects to delete
%         *regardless* of 'HandleVisibility'. Only these objects are deleted.
%         '!all' denotes all objects
%
% FORMAT spm_figure('Close',F)
% Closes figures (deletion without confirmation)
% Also closes the docking container if empty.
% F - 'Tag' string or figure number of figure to clear, defaults to gcf
% 
% FORMAT spm_figure('Print',F)
% F - [Optional] Figure to print. ('Tag' or figure number)
%     Defaults to figure 'Tag'ed as 'Graphics'.
%     If none found, uses CurrentFigure if available.
% If objects 'Tag'ed 'NextPage' and 'PrevPage' are found, then the
% pages are shown and printed in order. In brief, pages are held as
% separate axes, with ony one 'Visible' at any one time. The handles of
% the "page" axes are stored in the 'UserData' of the 'NextPage'
% object, while the 'PrevPage' object holds the current page number.
%
% FORMAT [hNextPage, hPrevPage, hPageNo] = spm_figure('NewPage',hPage)
% SPM pagination function: Makes objects with handles hPage paginated
% Creates pagination buttons if necessary.
% hPage                         - Handles of objects to stick to this page
% hNextPage, hPrevPage, hPageNo - Handles of pagination controls
%
% FORMAT spm_figure('TurnPage',move,F)
% SPM pagination function: Turn to specified page
%
% FORMAT spm_figure('DeletePageControls',F)
% SPM pagination function: Deletes page controls
% F - [Optional] Figure in which to attempt to turn the page
%         Defaults to 'Graphics' 'Tag'ged window
%
% FORMAT n = spm_figure('#page')
% Returns the current number of pages.
%
% FORMAT n = spm_figure('CurrentPage');
% Return the current page number.
%
% FORMAT spm_figure('WaterMark',F,str,Tag,Angle,Perm)
% Adds watermark to figure windows.
% F - Figure for watermark. Defaults to gcf
% str   - Watermark string. Defaults (missing or empty) to SPM
% Tag   - Tag for watermark axes. Defaults to ''
% Angle - Angle for watermark. Defaults to -45
% Perm  - If specified, then watermark is permanent (HandleVisibility 'off')
%
% FORMAT F = spm_figure('CreateWin',Tag,Name,Visible)
% Creates a full length WhiteBg figure 'Tag'ged Tag (if specified).
% F   - Figure created
% Tag     - Tag for window
% Name    - Name for window
% Visible - 'on' or 'off'
%
% FORMAT spm_figure('CreateBar',F)
% Creates toolbar in figure F (defaults to gcf). F can be a 'Tag'
%
% FORMAT spm_figure('ColorMap')
% Callback for "ColorMap" menu
%
% FORMAT spm_figure('FontSize')
% Callback for "FontSize" menu
%__________________________________________________________________________


%-Condition arguments
%--------------------------------------------------------------------------
if ~nargin, Action = 'Create'; else Action = varargin{1}; end

%==========================================================================
switch lower(Action), case 'create'
%==========================================================================
% F = spm_figure('Create',Tag,Name,Visible)

if nargin<4, Visible='on'; else Visible=varargin{4}; end
if nargin<3, Name=''; else Name=varargin{3}; end
if nargin<2, Tag=''; else Tag=varargin{2}; end

F = spm_figure('CreateWin',Tag,Name,Visible);
spm_figure('CreateBar',F);
%spm_figure('FigContextMenu',F);
varargout = {F};

%==========================================================================
case 'findwin'
%==========================================================================
% F=spm_figure('FindWin',F)
% F=spm_figure('FindWin',Tag)
%-Find window: Find window with FigureNumber# / 'Tag' attribute
%-Returns empty if window cannot be found - deletes multiple tagged figs.

if nargin<2, F='Graphics'; else F=varargin{2}; end

if isempty(F)
    % Leave F empty
elseif ischar(F)
    % Finds Graphics window with 'Tag' string - delete multiples
    Tag = F;
    F = findall(allchild(0),'Flat','Tag',Tag);
    if length(F) > 1
        % Multiple Graphics windows - close all but most recent
        close(F(2:end))
        F = F(1);
    end
else
    % F is supposed to be a figure number - check it
    if ~any(F==allchild(0)), F=[]; end
end
varargout = {F};

%==========================================================================
case 'getwin'
%==========================================================================
% F=spm_figure('GetWin',Tag)
%-Like spm_figure('FindWin',Tag), except that if no such 'Tag'ged figure
% is found and 'Tag' is recognized, one is created.

if nargin<2, Tag='Graphics'; else Tag=varargin{2}; end
F = spm_figure('FindWin',Tag);

if isempty(F)
    if ischar(Tag)
        switch Tag
            case 'Graphics'
                F = spm_figure('Create','Graphics','Graphics');
            case 'DEM'
                F = spm_figure('Create','DEM','Dynamic Expectation Maximisation');
            case 'DFP'
                F = spm_figure('Create','DFP','Variational filtering');
            case 'FMIN'
                F = spm_figure('Create','FMIN','Function minimisation');
            case 'MFM'
                F = spm_figure('Create','MFM','Mean-field and neural mass models');
            case 'MVB'
                F = spm_figure('Create','MVB','Multivariate Bayes');
            case 'SI'
                F = spm_figure('Create','SI','System Identification');
            case 'PPI'
                F = spm_figure('Create','PPI','Physio/Psycho-Physiologic Interaction');
            case 'Interactive'
                F = spm('CreateIntWin');
            case 'Satellite'
                F = spm_figure('CreateSatWin');
            otherwise
                F = spm_figure('Create',Tag,Tag);
        end
    end
else
    figure(F);
end
varargout = {F};

%==========================================================================
case 'select'
%==========================================================================
% spm_figure('Select',F)
F=varargin{2};
if ishandle(F)
    set(0, 'CurrentFigure', F)
end

%==========================================================================
case 'focus'
%==========================================================================
% spm_figure('Focus',F)
if nargin<2, F=get(0,'CurrentFigure'); else F=varargin{2}; end
if ishandle(F)
    figure(F);
end


%==========================================================================
case 'parentfig'
%==========================================================================
% F=spm_figure('ParentFig',h)

warning('spm_figure(''ParentFig'',h) is deprecated. Use ANCESTOR instead.');
if nargin<2, error('No object specified'), else h=varargin{2}; end
F = ancestor(h,'figure');
varargout = {F};

%==========================================================================
case 'clear'
%==========================================================================
% spm_figure('Clear',F,Tags)

%-Sort out arguments
if nargin<3, Tags=[]; else Tags=varargin{3}; end
if nargin<2, F=get(0,'CurrentFigure'); else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Clear figure
isdocked = strcmp(get(F,'WindowStyle'),'docked');
if isempty(Tags)
    %-Clear figure of objects with 'HandleVisibility' 'on'
    pos = get(F,'Position');
    delete(findall(allchild(F),'flat','HandleVisibility','on'));
    drawnow
    pause(0.05);
    if ~isdocked, set(F,'Position',pos); end
    %-Reset figures callback functions
    zoom(F,'off');
    pan(F,'off');
    rotate3d(F,'off');
    set(F,'KeyPressFcn','',...
        'WindowButtonDownFcn','',...
        'WindowButtonMotionFcn','',...
        'WindowButtonUpFcn','')
    %-If this is the 'Interactive' window, reset name & UserData
    if strcmp(get(F,'Tag'),'Interactive')
        set(F,'Name','','UserData',[]), end
else
    %-Clear specified objects from figure
    if ischar(Tags); Tags=cellstr(Tags); end
    if any(strcmp(Tags(:),'!all'))
        delete(allchild(F))
    else
        for tag = Tags(:)'
        delete(findall(allchild(F),'flat','Tag',tag{:}));
        end
    end 
end
set(F,'Pointer','Arrow')
%if ~isdocked && ~spm('CmdLine'), movegui(F); end

%==========================================================================
case 'close'
%==========================================================================
% spm_figure('Close',F)

%-Sort out arguments
if nargin < 2, F = gcf; else F = varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Detect if SPM windows are in docked mode
hMenu = spm_figure('FindWin','Menu');
isdocked = strcmp(get(hMenu,'WindowStyle'),'docked');

%-Close figures (and deleted without confirmation)
delete(F);

%-If in docked mode and closing SPM, close the container as well
if isdocked && ismember(hMenu,F)
    try
        desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
        group   = ['Statistical Parametric Mapping (' spm('Ver') ')'];
        hContainer = desktop.getGroupContainer(group);
        hContainer.getTopLevelAncestor.hide;
    end
end

%==========================================================================
case 'print'
%==========================================================================
% spm_figure('Print',F,fname,opt)

%-Arguments & defaults
if nargin<4, opt=spm_get_defaults('ui.print'); else opt=varargin{4}; end
if nargin<3, fname=''; else fname=varargin{3}; end
if nargin<2, F='Graphics'; else F=varargin{2}; end

%-Find window to print, default to gcf if specified figure not found
% Return if no figures
if ~isempty(F), F = spm_figure('FindWin',F); end
if  isempty(F), F = get(0,'CurrentFigure'); end
if  isempty(F), return, end

%-Note current figure, & switch to figure to print
cF = get(0,'CurrentFigure');
set(0,'CurrentFigure',F)

%-See if window has paging controls
hNextPage = findall(F,'Tag','NextPage');
hPrevPage = findall(F,'Tag','PrevPage');
hPageNo   = findall(F,'Tag','PageNo');
iPaged    = ~isempty(hNextPage);

%-Temporarily change all units to normalized prior to printing
H  = findall(allchild(F),'flat','Type','axes');
if ~isempty(H)
    un = cellstr(get(H,'Units'));
    set(H,'Units','normalized');
end

%-Print
if ~iPaged
    spm_print(fname,F,opt);
else
    hPg    = get(hNextPage,'UserData');
    Cpage  = get(hPageNo,  'UserData');
    nPages = size(hPg,1);

    set([hNextPage,hPrevPage,hPageNo],'Visible','off');
    if Cpage~=1
        set(hPg{Cpage,1},'Visible','off');
    end
    for p = 1:nPages
        set(hPg{p,1},'Visible','on');
        spm_print(fname,F,opt);
        set(hPg{p,1},'Visible','off');
    end
    set(hPg{Cpage,1},'Visible','on');
    set([hNextPage,hPrevPage,hPageNo],'Visible','on');
end
if ~isempty(H), set(H,{'Units'},un); end
set(0,'CurrentFigure',cF);

%==========================================================================
case 'printto'
%==========================================================================
%spm_figure('PrintTo',F)

%-Arguments & defaults
if nargin<2, F='Graphics'; else F=varargin{2}; end

%-Find window to print, default to gcf if specified figure not found
% Return if no figures
F=spm_figure('FindWin',F);
if isempty(F), F = get(0,'CurrentFigure'); end
if isempty(F), return, end

def = spm_get_defaults('ui.print');
if ischar(def), [def, i] = spm_print('format',def); end
pf  = spm_print('format');
pf  = pf([i setxor(1:numel(pf),i)]);
fs  = cell(numel(pf),2);
for i=1:numel(pf)
    fs{i,1} = ['*' pf(i).ext];
    fs{i,2} = [pf(i).name sprintf(' (*%s)',pf(i).ext)];
end
[fn, pn, fi] = uiputfile(fs,'Save figure as');
if isequal(fn,0) || isequal(pn,0), return, end

fname = fullfile(pn, fn);
spm_figure('Print',F,fname,pf(fi));

%==========================================================================
case 'newpage'
%==========================================================================
% [hNextPage, hPrevPage, hPageNo] = spm_figure('NewPage',h)

if nargin<2 || isempty(varargin{2}), error('No handles to paginate')
else h=varargin{2}(:)'; end

%-Work out which figure we're in
F = ancestor(h(1),'figure');

hNextPage = findall(F,'Tag','NextPage');
hPrevPage = findall(F,'Tag','PrevPage');
hPageNo   = findall(F,'Tag','PageNo');

%-Create pagination widgets if required
%--------------------------------------------------------------------------
if isempty(hNextPage)
    WS = spm('WinScale');
    FS = spm('FontSizes');
    SatFig = findall(0,'Tag','Satellite');
    if ~isempty(SatFig)
        SatFigPos    = get(SatFig,'Position');
        hNextPagePos = [SatFigPos(3)-25 15 15 15];
        hPrevPagePos = [SatFigPos(3)-40 15 15 15];
        hPageNo      = [SatFigPos(3)-40  5 30 10];
    else
        hNextPagePos = [580 022 015 015].*WS;
        hPrevPagePos = [565 022 015 015].*WS;
        hPageNo      = [550 005 060 015].*WS;
    end
    
    hNextPage = uicontrol(F,'Style','Pushbutton',...
        'HandleVisibility','on',...
        'String','>','FontSize',FS(10),...
        'ToolTipString','next page',...
        'Callback','spm_figure(''TurnPage'',''+1'',gcbf)',...
        'Position',hNextPagePos,...
        'ForegroundColor',[0 0 0],...
        'Tag','NextPage','UserData',[]);
    hPrevPage = uicontrol(F,'Style','Pushbutton',...
        'HandleVisibility','on',...
        'String','<','FontSize',FS(10),...
        'ToolTipString','previous page',...
        'Callback','spm_figure(''TurnPage'',''-1'',gcbf)',...
        'Position',hPrevPagePos,...
        'Visible','on',...
        'Enable','off',...
        'Tag','PrevPage');
    hPageNo = uicontrol(F,'Style','Text',...
        'HandleVisibility','on',...
        'String','1',...
        'FontSize',FS(6),...
        'HorizontalAlignment','center',...
        'BackgroundColor','w',...
        'Position',hPageNo,...
        'Visible','on',...
        'UserData',1,...
        'Tag','PageNo','UserData',1);
end

%-Add handles for this page to UserData of hNextPage
%-Make handles for this page invisible if PageNo>1
%--------------------------------------------------------------------------
mVis    = strcmp('on',get(h,'Visible'));
mHit    = strcmp('on',get(h,'HitTest'));
hPg     = get(hNextPage,'UserData');
if isempty(hPg)
    hPg = {h(mVis), h(~mVis), h(mHit), h(~mHit)};
else
    hPg = [hPg; {h(mVis), h(~mVis), h(mHit), h(~mHit)}];
    set(h(mVis),'Visible','off');
        set(h(mHit),'HitTest','off');
end
set(hNextPage,'UserData',hPg)

%-Return handles to pagination controls if requested
if nargout>0, varargout = {[hNextPage, hPrevPage, hPageNo]}; end

%==========================================================================
case 'turnpage'
%==========================================================================
% spm_figure('TurnPage',move,F)

if nargin<3, F='Graphics'; else F=varargin{3}; end
if nargin<2, move=1; else move=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findall(F,'Tag','NextPage');
hPrevPage = findall(F,'Tag','PrevPage');
hPageNo   = findall(F,'Tag','PageNo');
if isempty(hNextPage), return, end
hPg       = get(hNextPage,'UserData');
Cpage     = get(hPageNo,  'UserData');
nPages    = size(hPg,1);

%-Sort out new page number
if ischar(move), Npage = Cpage+eval(move); else Npage = move; end
Npage = max(min(Npage,nPages),1);

%-Make current page invisible, new page visible, set page number string
set(hPg{Cpage,1},'Visible','off');
set(hPg{Cpage,3},'HitTest','off');
set(hPg{Npage,1},'Visible','on');
set(hPg{Npage,3},'HitTest','on');
set(hPageNo,'UserData',Npage,'String',sprintf('%d / %d',Npage,nPages))

for k = 1:length(hPg{Npage,1})
    if strcmp(get(hPg{Npage,1}(k),'Type'),'axes')
        axes(hPg{Npage,1}(k));
    end
end

%-Disable appropriate page turning control if on first/last page
if Npage==1, set(hPrevPage,'Enable','off')
else set(hPrevPage,'Enable','on'), end
if Npage==nPages, set(hNextPage,'Enable','off')
else set(hNextPage,'Enable','on'), end

%==========================================================================
case 'deletepagecontrols'
%==========================================================================
% spm_figure('DeletePageControls',F)

if nargin<2, F='Graphics'; else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findall(F,'Tag','NextPage');
hPrevPage = findall(F,'Tag','PrevPage');
hPageNo   = findall(F,'Tag','PageNo');

delete([hNextPage hPrevPage hPageNo])

%==========================================================================
case '#page'
%==========================================================================
% n = spm_figure('#Page',F)

if nargin<2, F='Graphics'; else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findall(F,'Tag','NextPage');
if isempty(hNextPage)
    n = 1;
else
    n = size(get(hNextPage,'UserData'),1)+1;
end
varargout = {n};

%==========================================================================
case 'currentpage'
%==========================================================================
% n = spm_figure('CurrentPage', F)

if nargin<2, F='Graphics'; else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hPageNo   = findall(F,'Tag','PageNo');
Cpage     = get(hPageNo,  'UserData');

varargout = {Cpage};

%==========================================================================
case 'watermark'
%==========================================================================
% spm_figure('WaterMark',F,str,Tag,Angle,Perm)

if nargin<6, HVis='on'; else HVis='off'; end
if nargin<5, Angle=-45; else Angle=varargin{5}; end
if nargin<4 || isempty(varargin{4}), Tag = 'WaterMark'; else Tag=varargin{4}; end
if nargin<3 || isempty(varargin{3}), str = 'SPM';       else str=varargin{3}; end
if nargin<2, if any(allchild(0)), F=gcf; else F=''; end
else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Specify watermark color from background colour
Colour = get(F,'Color');
%-Only mess with grayscale backgrounds
if ~all(Colour==Colour(1)), return, end
%-Work out colour - lighter unless grey value > 0.9
Colour = Colour+(2*(Colour(1)<0.9)-1)*0.02;

cF = get(0,'CurrentFigure');
set(0,'CurrentFigure',F)
Units=get(F,'Units');
set(F,'Units','normalized');
h = axes('Position',[0.45,0.5,0.1,0.1],...
    'Units','normalized',...
    'Visible','off',...
    'Tag',Tag);
set(F,'Units',Units)
text(0.5,0.5,str,...
    'FontSize',spm('FontSize',80),...
    'FontWeight','Bold',...
    'FontName',spm_platform('Font','times'),...
    'Rotation',Angle,...
    'HorizontalAlignment','Center',...
    'VerticalAlignment','middle',...
    'Color',Colour,...
    'ButtonDownFcn',[...
        'if strcmp(get(gcbf,''SelectionType''),''open''),',...
            'delete(get(gcbo,''Parent'')),',...
        'end'])
set(h,'HandleVisibility',HVis)
set(0,'CurrentFigure',cF)

%==========================================================================
case 'createwin'
%==========================================================================
% F=spm_figure('CreateWin',Tag,Name,Visible)

if nargin<4 || isempty(varargin{4}), Visible='on'; else Visible=varargin{4}; end
if nargin<3, Name=''; else Name = varargin{3}; end
if nargin<2, Tag='';  else Tag  = varargin{2}; end

FS   = spm('FontSizes');          %-Scaled font sizes
PF   = spm_platform('fonts');     %-Font names (for this platform)
Rect = spm('WinSize','Graphics'); %-Graphics window rectangle
S0   = spm('WinSize','0',1);      %-Screen size (of the current monitor)

F    = figure(...
    'Tag',Tag,...
    'Position',[S0(1) S0(2) 0 0] + Rect,...
    'Resize','off',...
    'Color','w',...
    'ColorMap',gray(64),...
    'DefaultTextColor','k',...
    'DefaultTextInterpreter','none',...
    'DefaultTextFontName',PF.helvetica,...
    'DefaultTextFontSize',FS(10),...
    'DefaultAxesColor','w',...
    'DefaultAxesXColor','k',...
    'DefaultAxesYColor','k',...
    'DefaultAxesZColor','k',...
    'DefaultAxesFontName',PF.helvetica,...
    'DefaultPatchFaceColor','k',...
    'DefaultPatchEdgeColor','k',...
    'DefaultSurfaceEdgeColor','k',...
    'DefaultLineColor','k',...
    'DefaultUicontrolFontName',PF.helvetica,...
    'DefaultUicontrolFontSize',FS(10),...
    'DefaultUicontrolInterruptible','on',...
    'PaperType','A4',...
    'PaperUnits','normalized',...
    'PaperPosition',[.0726 .0644 .854 .870],...
    'InvertHardcopy','off',...
    'Renderer',spm_get_defaults('renderer'),...
    'Visible','off',...
    'Toolbar','none');
if ~isempty(Name)
    set(F,'Name',sprintf('%s: %s',spm('Version'),Name),'NumberTitle','off');
end
set(F,'Visible',Visible)
varargout = {F};

isdocked = strcmp(get(spm_figure('FindWin','Menu'),'WindowStyle'),'docked');
if isdocked
    try
        desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
        group   = ['Statistical Parametric Mapping (' spm('Ver') ')'];
        set(getJFrame(F),'GroupName',group);
        set(F,'WindowStyle','docked');
    end
end

%==========================================================================
 case 'createsatwin'
%==========================================================================
% F=spm_figure('CreateSatWin')
F = spm_figure('FindWin','Satellite');
if ~isempty(F)
    spm_figure('Focus', F)
else
    FS   = spm('FontSizes');             %-Scaled font sizes
    PF   = spm_platform('fonts');        %-Font names
    WS   = spm('WinSize','0','raw');     %-Screen size (of current monitor)
    Rect = [WS(1)+5 WS(4)*.40 WS(3)*.49 WS(4)*.57];
    F = figure(...
        'Tag','Satellite',...
        'Position',Rect,...
        'Resize','off',...
        'MenuBar','none',...
        'Name','SPM: Satellite Results Table',...
        'Numbertitle','off',...
        'Color','w',...
        'ColorMap',gray(64),...
        'DefaultTextColor','k',...
        'DefaultTextInterpreter','none',...
        'DefaultTextFontName',PF.helvetica,...
        'DefaultTextFontSize',FS(10),...
        'DefaultAxesColor','w',...
        'DefaultAxesXColor','k',...
        'DefaultAxesYColor','k',...
        'DefaultAxesZColor','k',...
        'DefaultAxesFontName',PF.helvetica,...
        'DefaultPatchFaceColor','k',...
        'DefaultPatchEdgeColor','k',...
        'DefaultSurfaceEdgeColor','k',...
        'DefaultLineColor','k',...
        'DefaultUicontrolFontName',PF.helvetica,...
        'DefaultUicontrolFontSize',FS(10),...
        'DefaultUicontrolInterruptible','on',...
        'PaperType','A4',...
        'PaperUnits','normalized',...
        'PaperPosition',[.0726 .0644 .854 .870],...
        'InvertHardcopy','off',...
        'Renderer',spm_get_defaults('renderer'),...
        'Visible','on');
end
varargout = {F};

%==========================================================================
 case 'createbar'
%==========================================================================
% spm_figure('CreateBar',F)

if nargin<2, if any(allchild(0)), F=gcf; else F=''; end
else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Help Menu
t0 = findall(allchild(F),'Flat','Label','&Help'); 
delete(allchild(t0)); set(t0,'Callback',''); set(t0,'Tag','');
if isempty(t0), t0 = uimenu( F,'Label','&Help'); end
pos = get(t0,'Position');
uimenu(t0,'Label','SPM Help','CallBack','spm_help');
uimenu(t0,'Label','SPM Manual (PDF)',...
    'CallBack','try,open(fullfile(spm(''dir''),''man'',''manual.pdf''));end');
t1=uimenu(t0,'Label','SPM &Web Resources');
uimenu(t1,'Label','SPM Web &Site',...
    'CallBack','web(''http://www.fil.ion.ucl.ac.uk/spm/'');');
uimenu(t1,'Label','SPM &WikiBook',...
    'CallBack','web(''http://en.wikibooks.org/wiki/SPM'');');
uimenu(t1,'Separator','on','Label','SPM &Extensions',...
    'CallBack','web(''http://www.fil.ion.ucl.ac.uk/spm/ext/'');');

%-Check Menu
if ~isdeployed
    uimenu(t0,'Separator','on','Label','SPM Check Installation',...
        'CallBack','spm_check_installation(''full'')');
    uimenu(t0,'Label','SPM Check for Updates',...
        'CallBack',@spm_check_update);
end

%- About Menu
uimenu(t0,'Separator','on','Label',['&About ' spm('Ver')],...
    'CallBack',@spm_about);
if strcmpi(spm_check_version,'matlab')
    uimenu(t0,'Label','&About MATLAB',...
        'CallBack','web(''http://www.mathworks.com/matlab/'');');
else
    uimenu(t0,'Label','&About GNU Octave',...
        'CallBack','web(''http://www.octave.org/'');');
end

%-Figure Menu
t0=uimenu(F, 'Position',pos, 'Label','&SPM Figure', 'HandleVisibility','off', 'Callback',@myfigmenu);

%-Show All Figures
uimenu(t0, 'Label','Show All &Windows', 'HandleVisibility','off',...
    'CallBack','spm(''Show'');');

if strcmpi(spm_check_version,'matlab')
    %-Show MATLAB Command Window
    uimenu(t0, 'Label','Show &MATLAB Window', 'HandleVisibility','off',...
        'CallBack','commandwindow;');

    %-Dock SPM Figures
    uimenu(t0, 'Label','&Dock SPM Windows', 'HandleVisibility','off',...
        'CallBack',@mydockspm);
end

%-Print Menu
%t1=uimenu(t0, 'Label','&Save Figure', 'HandleVisibility','off','Separator','on');
uimenu(t0,    'Label','&Save Figure', 'HandleVisibility','off', ...
    'CallBack','spm_figure(''Print'',gcbf)', 'Separator','on');
uimenu(t0,    'Label','Save Figure &As...', 'HandleVisibility','off', ...
    'CallBack','spm_figure(''PrintTo'',gcbf)');

%-Copy Figure
if ispc
    uimenu(t0, 'Label','Co&py Figure', 'HandleVisibility','off',...
       'CallBack','editmenufcn(gcbf,''EditCopyFigure'')');
end

%-Clear Menu
uimenu(t0,    'Label','&Clear Figure', 'HandleVisibility','off', ...
    'CallBack','spm_figure(''Clear'',gcbf)');

%-Close non-SPM figures
uimenu(t0,    'Label','C&lose non-SPM Figures', 'HandleVisibility','off', ...
    'CallBack',@myclosefig);

%-Colour Menu
t1=uimenu(t0, 'Label','C&olours',  'HandleVisibility','off','Separator','on');
t2=uimenu(t1, 'Label','Colormap');
uimenu(t2,    'Label','Gray',      'CallBack','spm_figure(''ColorMap'',''gray'')');
uimenu(t2,    'Label','Hot',       'CallBack','spm_figure(''ColorMap'',''hot'')');
uimenu(t2,    'Label','Pink',      'CallBack','spm_figure(''ColorMap'',''pink'')');
uimenu(t2,    'Label','Jet',       'CallBack','spm_figure(''ColorMap'',''jet'')');
uimenu(t2,    'Label','Gray-Hot',  'CallBack','spm_figure(''ColorMap'',''gray-hot'')');
uimenu(t2,    'Label','Gray-Cool', 'CallBack','spm_figure(''ColorMap'',''gray-cool'')');
uimenu(t2,    'Label','Gray-Pink', 'CallBack','spm_figure(''ColorMap'',''gray-pink'')');
uimenu(t2,    'Label','Gray-Jet',  'CallBack','spm_figure(''ColorMap'',''gray-jet'')');
t2=uimenu(t1, 'Label','Effects');
uimenu(t2,    'Label','Invert',    'CallBack','spm_figure(''ColorMap'',''invert'')');
uimenu(t2,    'Label','Brighten',  'CallBack','spm_figure(''ColorMap'',''brighten'')');
uimenu(t2,    'Label','Darken',    'CallBack','spm_figure(''ColorMap'',''darken'')');

%-Font Size Menu
t1=uimenu(t0, 'Label','&Font Size', 'HandleVisibility','off');
uimenu(t1, 'Label','&Increase', 'CallBack','spm_figure(''FontSize'',1)',  'Accelerator', '=');
uimenu(t1, 'Label','&Decrease', 'CallBack','spm_figure(''FontSize'',-1)', 'Accelerator', '-');

%-Renderer Menu
t1=uimenu(t0, 'Label','R&enderer', 'HandleVisibility','off');
uimenu(t1, 'Label', 'painters', 'CallBack','spm_get_defaults(''renderer'',''painters'');set(gcf,''Renderer'',''painters'');');
uimenu(t1, 'Label', 'zbuffer',  'CallBack','spm_get_defaults(''renderer'',''zbuffer'');set(gcf,''Renderer'',''zbuffer'');');
uimenu(t1, 'Label', 'OpenGL',   'CallBack','spm_get_defaults(''renderer'',''opengl'');set(gcf,''Renderer'',''opengl'');');

%-Satellite Table
uimenu(t0,    'Label','&Results Table', 'HandleVisibility','off', ...
    'Separator','on', 'Callback','spm_figure(''CreateSatWin'');');
    
%==========================================================================
case 'figcontextmenu'
%==========================================================================
% h = spm_figure('FigContextMenu',F)

if nargin<2
    F = get(0,'CurrentFigure');
    if isempty(F), error('no figure'), end
else
    F = spm_figure('FindWin',varargin{2});
    if isempty(F), error('no such figure'), end
end
h     = uicontextmenu('Parent',F,'HandleVisibility','CallBack');
copy_menu(F,h);
set(F,'UIContextMenu',h)
varargout = {h};

%==========================================================================
case 'colormap'
%==========================================================================
% spm_figure('ColorMap',ColAction)

if nargin<2, ColAction='gray'; else ColAction=varargin{2}; end

switch lower(ColAction), case 'gray'
    colormap(gray(64))
case 'hot'
    colormap(hot(64))
case 'pink'
    colormap(pink(64))
case 'jet'
    colormap(jet(64))
case 'gray-hot'
    tmp = hot(64 + 16);  tmp = tmp((1:64) + 16,:);
    colormap([gray(64); tmp]);
case 'gray-cool'
    cool = [zeros(10,1) zeros(10,1) linspace(0.5,1,10)';
            zeros(31,1) linspace(0,1,31)' ones(31,1);
            linspace(0,1,23)' ones(23,1) ones(23,1) ];
    colormap([gray(64); cool]);
case 'gray-pink'
    tmp = pink(64 + 16); tmp = tmp((1:64) + 16,:);
    colormap([gray(64); tmp]);
case 'gray-jet'
    colormap([gray(64); jet(64)]);
case 'invert'
    colormap(flipud(colormap));
case 'brighten'
    colormap(brighten(colormap, 0.2));
case 'darken'
    colormap(brighten(colormap, -0.2));
otherwise
    error('Illegal ColAction specification');
end

%==========================================================================
case 'fontsize'
%==========================================================================
% spm_figure('FontSize',sz)

if nargin<2, sz=0; else sz=varargin{2}; end

h  = [get(0,'CurrentFigure') spm_figure('FindWin','Satellite')];
h  = [findall(h,'type','text'); findall(h,'type','uicontrol')];
fs = get(h,'fontsize');
if ~iscell(fs), fs = {fs}; end
if ~isempty(fs)
    set(h,{'fontsize'},cellfun(@(x) max(x+sz,eps),fs,'UniformOutput',false));
end


%==========================================================================
otherwise
%==========================================================================
warning(['Illegal Action string: ',Action])
end
return;


%==========================================================================
function myfigmenu(obj,evt)
%==========================================================================
hr = findall(obj,'Label','&Results Table');
try
    evalin('base','xSPM;');
    set(hr,'Enable','on');
catch
    set(hr,'Enable','off');
end
SatWindow = spm_figure('FindWin','Satellite');
if ~isempty(SatWindow)
    set(hr,'Checked','on');
else
    set(hr,'Checked','off');
end
rend = get(ancestor(obj,'figure'),'Renderer');
hr   = get(findall(obj,'Label','R&enderer'),'Children');
set(hr,'Checked','off');
set(hr(ismember(lower(get(hr,'Label')),rend)),'Checked','on');

% for compatibility when opening a FIG-file from a previous version of SPM
function myisresults(obj,evt)
return;
function mysatfig(obj,evt)
return;


%==========================================================================
function mydockspm(obj,evt)
%==========================================================================
% Largely inspired by setFigDockGroup from Yair Altman
% http://www.mathworks.com/matlabcentral/fileexchange/16650
hMenu = spm_figure('FindWin','Menu');
hInt  = spm_figure('FindWin','Interactive');
hGra  = spm_figure('FindWin','Graphics');
h     = [hMenu hInt hGra];
group   = ['Statistical Parametric Mapping (' spm('Ver') ')'];
try
    desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    if ~ismember(group,cell(desktop.getGroupTitles))
        desktop.addGroup(group);
    end
    for i=1:length(h)
        set(getJFrame(h(i)),'GroupName',group);
    end
    hContainer = desktop.getGroupContainer(group);
    set(hContainer,'userdata',group);
end
set(h,'WindowStyle','docked');
try, pause(0.5), desktop.setGroupDocked(group,false); end

%==========================================================================
function myclosefig(obj,evt)
%==========================================================================
hMenu = spm_figure('FindWin','Menu');
hInt  = spm_figure('FindWin','Interactive');
hGra  = spm_figure('FindWin','Graphics');
hSat  = spm_figure('FindWin','Satellite');
hBat  = spm_figure('FindWin','cfg_ui');
h     = setdiff(findobj(get(0,'children'),'flat','visible','on'), ...
    [hMenu; hInt; hGra; hSat; hBat; gcf]);
close(h,'force');

%==========================================================================
function copy_menu(F,G)
%==========================================================================
handles = findall(allchild(F),'Flat','Type','uimenu','Visible','on');
if isempty(handles), return; end;
for F1=handles(:)'
    if ~ismember(get(F1,'Label'),{'&Window' '&Desktop'})
        G1 = uimenu(G,'Label',get(F1,'Label'),...
            'CallBack',get(F1,'CallBack'),...
            'Position',get(F1,'Position'),...
            'Separator',get(F1,'Separator'));
        copy_menu(F1,G1);
    end
end

%==========================================================================
function jframe = getJFrame(h)
%==========================================================================
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
hhFig    = handle(h);
jframe   = [];
maxTries = 16;
while maxTries > 0
    try
        jframe = get(hhFig,'javaframe');
        if ~isempty(jframe)
            break;
        else
            maxTries = maxTries - 1;
            drawnow; pause(0.1);
        end
    catch
        maxTries = maxTries - 1;
        drawnow; pause(0.1);
    end
end
if isempty(jframe)
    error('Cannot retrieve the java frame from handle.');
end

%==========================================================================
function spm_about(obj,evt)
%==========================================================================
h = figure('MenuBar','none',...
           'NumberTitle','off',...
           'Name',['About ' spm('Version')],...
           'Resize','off',...
           'Toolbar','none',...
           'Tag','AboutSPM',...
           'WindowStyle','Modal',...
           'Color',[1 1 1],...
           'Visible','off',...
           'DoubleBuffer','on');
pos = get(h,'Position');
pos([3 4]) = [300 400];
set(h,'Position',pos);
set(h,'Visible','on');

a = axes('Parent',h, 'Units','pixels', 'Position',[50 201 200 200],...
    'Visible','off');
IMG = imread(fullfile(spm('Dir'),'help','images','spm12.png'));
image(IMG,'Parent',a); set(a,'Visible','off');

a = axes('Parent',h,'Units','pixels','Position',[0 0 300 400],...
    'Visible','off','Tag','textcont');
text(0.5,0.45,'Statistical Parametric Mapping','Parent',a,...
    'HorizontalAlignment','center','Color',[0 0 0],'FontWeight','Bold');
text(0.5,0.40,spm('Version'),'Parent',a,'HorizontalAlignment','center',...
    'Color',[1 1 1]);
text(0.5,0.30,'Wellcome Trust Centre for Neuroimaging','Parent',a,...
    'HorizontalAlignment','center','Color',[0 0 0],'FontWeight','Bold');
text(0.5,0.25,['Copyright (C) 1991,1994-' datestr(now,'yyyy')],...
    'Parent',a,'HorizontalAlignment','center','Color',[0 0 0]);
text(0.5,0.20,'http://www.fil.ion.ucl.ac.uk/spm/','Parent',a,...
    'HorizontalAlignment','center','Color',[0 0 0],...
    'ButtonDownFcn','web(''http://www.fil.ion.ucl.ac.uk/spm/'');');
if isdeployed
    text(0.5,0.15,['MATLAB(r). (c) 1984-' datestr(now,'yyyy') ' The MathWorks, Inc.'],...
    'Parent',a,'HorizontalAlignment','center','Color',[0 0 0]);
end

uicontrol('Style','pushbutton','String','Credits','Position',[40 25 60 25],...
    'Callback',@myscroll,'BusyAction','Cancel');
uicontrol('Style','pushbutton','String','OK','Position',[200 25 60 25],...
    'Callback','close(gcf)','BusyAction','Cancel');

%==========================================================================
function myscroll(obj,evt)
%==========================================================================
ax = findobj(gcf,'Tag','textcont');
cla(ax);
[current, previous] = spm_authors;
authors = {['*' spm('Ver') ' Developers*'] current{:} '' ...
           '*Collaborators and*' '*former Contributors*' previous{:} '' ...
           '*Thanks to the SPM community*'};
x = 0.2;
h = [];
for i=1:numel(authors)
    h(i) = text(0.5,x,authors{i},'Parent',ax,...
        'HorizontalAlignment','center','Color',col(x));
    if any(authors{i} == '*')
        set(h(i),'String',strrep(authors{i},'*',''),'FontWeight','Bold');
    end
    x = x - 0.05;
end
pause(0.5);
try
    for j=1:fix((0.5-(0.2-numel(authors)*0.05))/0.01)
        for i=1:numel(h)
            p  = get(h(i),'Position');
            p2 = p(2)+0.01;
            set(h(i),'Position',[p(1) p2 p(3)],'Color',col(p2));
            if p2 > 0.5, set(h(i),'Visible','off'); end
        end
        pause(0.1)
    end
end

%==========================================================================
function c = col(x)
%==========================================================================
if x < 0.4 && x > 0.3
    c = [1 1 1];
elseif x <= 0.3
    c = [1 1 1] - 6*abs(0.3-x);
else
    c = [1 1 1] - 6*abs(0.4-x);
end
c(c<0) = 0; c(c>1) = 1;
c = 1 - c;

%==========================================================================
function spm_check_update(obj,evt)
%==========================================================================
[sts, msg] = spm_update;
if ~isfinite(sts) || sts == 0
    spm('alert"',msg,'Update');
    return
end

str = {'SPM updates are available!','Do you want to install them now?'};
if spm_input(str,1,'bd','Yes|No',[1,0],1,'Update')
    spm_update(true);
else
    fprintf(msg);
end
