function varargout=spm_figure(varargin)
% Setup and callback functions for Graphics window
% FORMAT varargout=spm_figure(varargin)
%       - An embedded callback, multi-function function
%       - For detailed programmers comments, see format specifications
%         in main body of code
%_______________________________________________________________________
%
% spm_figure creates and manages the 'Graphics' window. This window and
% these facilities may be used independently of SPM, and any number of
% Graphics windows my be used within the same MATLAB session. (Though
% only one SPM 'Graphics' 'Tag'ed window is permitted.
%
% The Graphics window is provided with a menu bar at the top that
% facilitates editing and printing of the current graphic display,
% enabling interactive editing of graphic output prior to printing
% (e.g. selection of color maps, deleting, moving and editing graphics
% objects or adding text). (This menu is also provided as a figure
% background "ContextMenu" - right-clicking on the figure background
% should bring up the menu.)
%
% Print: Graphics windows with multi-page axes are printed page by page.
%
% Clear: Clears the Graphics window. If in SPM usage (figure 'Tag'ed as
% 'Graphics') then all SPM windows are cleared and reset.
%
% Colormap options:
% * gray, hot, pink: Sets the colormap to its default values and loads
%                 either a grayscale, 'hot metal' or color map.
% * gray-hot, etc: Creates a 'split' colormap {128 x 3 matrix}.
%          The lower half is a gray scale and the upper half
%          is 'hot metal' or 'pink'.  This color map is used for
%      viewing 'rendered' SPMs on a PET, MRI or other background images
%
% Colormap effects:
% * Invert: Inverts (flips) the current color map.
% * Brighten and Darken: Brighten and Darken the current colourmap
%      using the MATLAB BRIGHTEN command, with  beta's of +0.2 and -0.2
%      respectively.
%
% Editing: Right button ('alt' button) cancels operations
% * Cut  : Deletes the graphics object next selected (if deletable)
%          Select with middle mouse button to delete blocks of text,
%          or to delete individual elements from a plot.
% * Move : To re-position a text, uicontrol or axis object using a
%          'drag and drop' implementation (i.e. depress - move - release)
%          Using the middle 'extend' mouse button on a text object moves
%          the axes containing the text - i.e. blocks of text.
% * Size : Re-sizes the text, uicontrol or axis object next selected
%          {left button - decrease, middle button  - increase} by a factor
%          of 1.24 (or increase/decrease FontSize by 2 dpi)
% * Text : Creates an editable text widget that produces a text object as
%          its CallBack.
%          The text object is provided with a ContextMenu, obtained by
%          right-clicking ('alt') on the text, allowing text attributes
%          to be changed. Alternatively, the edit facilities on the window
%          menu bar or ContextMenu can be used.
% * Edit : To edit text, select a text object with the circle cursor,
%          and edit the text in the editable text widget that appears.
%          A middle 'extend' mouse click places a context menu on the text
%          object, facilitating easy modification of text atributes.
%
% For SPM usage, the figure should be 'Tag'ed as 'Graphics'.
%
% For SPM power users, and programmers, spm_figure provides utility
% routines for using the SPM graphics interface. Of particular use are
% the GetWin, FindWin and Clear functions See the embedded callback
% reference in the main body of spm_figure, below the help text.
%
% See also: spm_print, spm_clf
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id$


%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take. Recall   )
%( MATLAB's command-function duality: `spm_figure Create` is           )
%( equivalent to `spm_figure('Create')`.                               )
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
% FORMAT F = spm_figure('ParentFig',h)
% Finds window containing the object whose handle is specified
% h - Handle of object whose parent figure is required
%   - If a vector, then first object handle is used
% F - Number or parent figure
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
%
% FORMAT spm_figure('Print',F)
% F - [Optional] Figure to print. ('Tag' or figure number)
%     Defaults to figure 'Tag'ed as 'Graphics'.
%     If none found, uses CurrentFigure if avaliable.
% If objects 'Tag'ed 'NextPage' and 'PrevPage' are found, then the
% pages are shown and printed in order. In breif, pages are held as
% seperate axes, with ony one 'Visible' at any one time. The handles of
% the "page" axes are stored in the 'UserData' of the 'NextPage'
% object, while the 'PrevPage' object holds the current page number.
% See spm_help('!Disp') for details on setting up paging axes.
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
% FORMAT WS = spm_figure('GetWinScale')
% Returns ratios of current display dimensions to that of a 1152 x 900
% Sun display. WS=[Xratio,Yratio,Xratio,Yratio]. Used for scaling other
% GUI elements.
% (Function duplicated in spm.m, repeated to reduce inter-dependencies.)
%
% FORMAT FS = spm_figure('FontSizes',FS)
% Returns fontsizes FS scaled for the current display.
% FS     - (vector of) Font sizes to scale
%          [default [08,09,11,13,14,6:36]]
%
% FORMAT spm_figure('CreateBar',F)
% Creates toolbar in figure F (defaults to gcf). F can be a 'Tag'
% If the figure is 'Tag'ed as 'Graphics' (SPM usage), then the Print button
% callback is set to attempt to clear an 'Interactive' figure too.
%
% FORMAT spm_figure('ColorMap')
% Callback for "ColorMap" buttons
%
% FORMAT h = spm_figure('GraphicsHandle',F)
% GUI choose object for handle identification. LeftMouse 'normal' returns
% handle, MiddleMouse 'extend' returns parents handle, RightMouse 'alt' cancels.
% F - figure to do a GUI "handle ID" in [Default gcbf]
%_______________________________________________________________________


%-Condition arguments
%-----------------------------------------------------------------------
if (nargin==0), Action = 'Create'; else Action = varargin{1}; end

switch lower(Action), case 'create'
%=======================================================================
% F = spm_figure('Create',Tag,Name,Visible)
%-Condition arguments
if nargin<4, Visible='on'; else Visible=varargin{4}; end
if nargin<3, Name=''; else, Name=varargin{3}; end
if nargin<2, Tag=''; else Tag=varargin{2}; end

F = spm_figure('CreateWin',Tag,Name,Visible);
spm_figure('CreateBar',F);
spm_figure('FigContextMenu',F);
varargout = {F};


case 'findwin'
%=======================================================================
% F=spm_figure('FindWin',F)
% F=spm_figure('FindWin',Tag)
%-Find window: Find window with FigureNumber# / 'Tag' attribute
%-Returns empty if window cannot be found - deletes multiple tagged figs.

if nargin<2, F='Graphics'; else F=varargin{2}; end

shh = get(0,'showhiddenhandles');
set(0,'showhiddenhandles','on');

if isempty(F)
    % Leave F empty
elseif ischar(F)
    % Finds Graphics window with 'Tag' string - delete multiples
    Tag=F;
    F = findobj(get(0,'Children'),'Flat','Tag',Tag);
    if length(F) > 1
        % Multiple Graphics windows - close all but most recent
        close(F(2:end))
        F = F(1);
    end
else
    % F is supposed to be a figure number - check it
    if ~any(F==get(0,'Children')), F=[]; end
end
set(0,'showhiddenhandles',shh);
varargout = {F};


case 'getwin'
%=======================================================================
% F=spm_figure('GetWin',Tag)

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
            otherwise
                F = spm_figure('Create',Tag,Tag);
        end
    end
else
    set(0,'CurrentFigure',F);
    figure(F);
end
varargout = {F};

case 'parentfig'
%=======================================================================
% F=spm_figure('ParentFig',h)
if nargin<2, error('No object specified'), else h=varargin{2}; end
F = get(h(1),'Parent');
while ~strcmp(get(F,'Type'),'figure'), F=get(F,'Parent'); end
varargout = {F};


case 'clear'
%=======================================================================
% spm_figure('Clear',F,Tags)

%-Sort out arguments
%-----------------------------------------------------------------------
if nargin<3, Tags=[]; else Tags=varargin{3}; end
if nargin<2, F=get(0,'CurrentFigure'); else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Clear figure
%-----------------------------------------------------------------------
if isempty(Tags)
    %-Clear figure of objects with 'HandleVisibility' 'on'
    pos = get(F,'Position');
    delete(findobj(get(F,'Children'),'flat','HandleVisibility','on'));
    drawnow
    set(F,'Position',pos);

    %-Reset figures callback functions
    zoom(F,'off');
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
    cSHH = get(0,'ShowHiddenHandles');
    set(0,'ShowHiddenHandles','on')
    if ischar(Tags); Tags=cellstr(Tags); end
    if any(strcmp(Tags(:),'!all'))
        delete(get(F,'Children'))
    else
        for tag = Tags(:)'
        delete(findobj(get(F,'Children'),'flat','Tag',tag{:}));
        end
    end 
    set(0,'ShowHiddenHandles',cSHH)
end
set(F,'Pointer','Arrow')
movegui(F);


case 'print'
%=======================================================================
% spm_figure('Print',F,fname)

%-Arguments & defaults
if nargin<3, fname=''; else fname=varargin{3};end
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
hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');
iPaged    = ~isempty(hNextPage);

%-Construct print command
%-----------------------------------------------------------------------

%-Temporarily change all units to normalized prior to printing
% (Fixes bizzarre problem with stuff jumping around!)
%-----------------------------------------------------------------------
H  = findobj(get(F,'Children'),'flat','Type','axes');
if ~isempty(H),
    un = cellstr(get(H,'Units'));
    set(H,'Units','normalized')
end;

%-Print
%-----------------------------------------------------------------------
if ~iPaged
    spm_print(fname)
else
    hPg       = get(hNextPage,'UserData');
    Cpage     = get(hPageNo,  'UserData');
    nPages    = size(hPg,1);

    set([hNextPage,hPrevPage,hPageNo],'Visible','off')
    if Cpage~=1
        set(hPg{Cpage,1},'Visible','off'), end
    for p = 1:nPages
        set(hPg{p,1},'Visible','on');
        spm_print(fname);
        set(hPg{p,1},'Visible','off')
    end
    set(hPg{Cpage,1},'Visible','on')
    set([hNextPage,hPrevPage,hPageNo],'Visible','on')
end
if ~isempty(H), set(H,{'Units'},un); end;
set(0,'CurrentFigure',cF)


case 'printto'
%=======================================================================
%spm_figure('PrintTo',F)

%-Arguments & defaults
if nargin<2, F='Graphics'; else F=varargin{2}; end

%-Find window to print, default to gcf if specified figure not found
% Return if no figures
F=spm_figure('FindWin',F);
if isempty(F), F = get(0,'CurrentFigure'); end
if isempty(F), return, end

[fn pn] = uiputfile('*.ps', 'Print to File');
if fn == 0
    return;
end;

psname = fullfile(pn, fn);
spm_figure('Print',F,psname);

case 'newpage'
%=======================================================================
% [hNextPage, hPrevPage, hPageNo] = spm_figure('NewPage',h)
if nargin<2 || isempty(varargin{2}) error('No handles to paginate')
else h=varargin{2}(:)'; end

%-Work out which figure we're in
F = spm_figure('ParentFig',h(1));

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');

%-Create pagination widgets if required
%-----------------------------------------------------------------------
if isempty(hNextPage)
    WS = spm('WinScale');
    FS = spm('FontSizes');
    SatFig = findobj('Tag','Satellite');
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
%-----------------------------------------------------------------------
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


case 'turnpage'
%=======================================================================
% spm_figure('TurnPage',move,F)
if nargin<3, F='Graphics'; else F=varargin{3}; end
if nargin<2, move=1; else move=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');
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

for k = 1:length(hPg{Npage,1}) % VG
    if strcmp(get(hPg{Npage,1}(k),'Type'),'axes')
        axes(hPg{Npage,1}(k));
    end;
end;

%-Disable appropriate page turning control if on first/last page (for neatness)
if Npage==1, set(hPrevPage,'Enable','off')
else set(hPrevPage,'Enable','on'), end
if Npage==nPages, set(hNextPage,'Enable','off')
else set(hNextPage,'Enable','on'), end



case 'deletepagecontrols'
%=======================================================================
% spm_figure('DeletePageControls',F)
if nargin<2, F='Graphics'; else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
hPrevPage = findobj(F,'Tag','PrevPage');
hPageNo   = findobj(F,'Tag','PageNo');

delete([hNextPage hPrevPage hPageNo])


case '#page'
%=======================================================================
% n = spm_figure('#Page',F)
if nargin<2, F='Graphics'; else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hNextPage = findobj(F,'Tag','NextPage');
if isempty(hNextPage)
    n = 1;
else
    n = size(get(hNextPage,'UserData'),1)+1;
end
varargout = {n};

case 'currentpage'
%=======================================================================
% n = spm_figure('CurrentPage', F)
if nargin<2, F='Graphics'; else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('No Graphics window'), end

hPageNo   = findobj(F,'Tag','PageNo');
Cpage     = get(hPageNo,  'UserData');

varargout = {Cpage};

case 'watermark'
%=======================================================================
% spm_figure('WaterMark',F,str,Tag,Angle,Perm)
if nargin<6, HVis='on'; else HVis='off'; end
if nargin<5, Angle=-45; else Angle=varargin{5}; end
if nargin<4 || isempty(varargin{4}), Tag = 'WaterMark'; else Tag=varargin{4}; end
if nargin<3 || isempty(varargin{3}), str = 'SPM';       else str=varargin{3}; end
if nargin<2, if any(get(0,'Children')), F=gcf; else F=''; end
else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Specify watermark color from background colour
%-----------------------------------------------------------------------
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
    'EraseMode','normal',...
    'Color',Colour,...
    'ButtonDownFcn',[...
        'if strcmp(get(gcbf,''SelectionType''),''open''),',...
            'delete(get(gcbo,''Parent'')),',...
        'end'])
set(h,'HandleVisibility',HVis)
set(0,'CurrentFigure',cF)


case 'createwin'
%=======================================================================
% F=spm_figure('CreateWin',Tag,Name,Visible)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<4 || isempty(varargin{4}), Visible='on'; else Visible=varargin{4}; end
if nargin<3, Name=''; else Name = varargin{3}; end
if nargin<2, Tag='';  else Tag  = varargin{2}; end

WS   = spm('WinScale');           %-Window scaling factors
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
    set(F,'Name',sprintf('%s%s: %s',spm('ver'),...
        spm('GetUser',' (%s)'),Name),'NumberTitle','off')
end
set(F,'Visible',Visible)
varargout = {F};


case 'getwinscale'
%=======================================================================
% WS = spm_figure('GetWinScale')
warning('spm_figure(''GetWinScale''... is Grandfathered: use spm(''WinScale''')
varargout = {spm('WinScale')};


case 'fontsizes'
%=======================================================================
% FS = spm_figure('FontSizes',FS)
warning('spm_figure(''FontSizes''... is Grandfathered: use spm(''FontSizes''')
if nargin<2, FS=[08,09,11,13,14,6:36]; else FS=varargin{2}; end
varargout = {round(FS*min(spm('WinScale')))};


%=======================================================================
 case 'createbar'
%=======================================================================
% spm_figure('CreateBar',F)
if nargin<2, if any(get(0,'Children')), F=gcf; else F=''; end
else F=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

cSHH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

t0 = findobj(get(F,'Children'),'Flat','Label','&Help');
if isempty(t0), t0 = uimenu( F,'Label','&Help'); end;

set(findobj(t0,'Position',1),'Separator','on');
uimenu(t0,'Position',1,...
    'Label','SPM web',...
    'CallBack','web(''http://www.fil.ion.ucl.ac.uk/spm/'');');
uimenu(t0,'Position',1,...
    'Label','SPM help','ForegroundColor',[0 1 0],...
    'CallBack','spm_help');

t0=uimenu( F,'Label','Colours','HandleVisibility','off');
t1=uimenu(t0,'Label','ColorMap');
uimenu(t1,'Label','Gray',     'CallBack','spm_figure(''ColorMap'',''gray'')');
uimenu(t1,'Label','Hot',      'CallBack','spm_figure(''ColorMap'',''hot'')');
uimenu(t1,'Label','Pink',     'CallBack','spm_figure(''ColorMap'',''pink'')');
uimenu(t1,'Label','Jet','CallBack','spm_figure(''ColorMap'',''jet'')');
uimenu(t1,'Label','Gray-Hot', 'CallBack','spm_figure(''ColorMap'',''gray-hot'')');
uimenu(t1,'Label','Gray-Cool','CallBack','spm_figure(''ColorMap'',''gray-cool'')');
uimenu(t1,'Label','Gray-Pink','CallBack','spm_figure(''ColorMap'',''gray-pink'')');
uimenu(t1,'Label','Gray-Jet', 'CallBack','spm_figure(''ColorMap'',''gray-jet'')');
t1=uimenu(t0,'Label','Effects');
uimenu(t1,'Label','Invert','CallBack','spm_figure(''ColorMap'',''invert'')');
uimenu(t1,'Label','Brighten','CallBack','spm_figure(''ColorMap'',''brighten'')');
uimenu(t1,'Label','Darken','CallBack','spm_figure(''ColorMap'',''darken'')');
uimenu( F,'Label','Clear','HandleVisibility','off','CallBack','spm_figure(''Clear'',gcbf)');
t0=uimenu( F,'Label','SPM-Print','HandleVisibility','off');
%uimenu( F,'Label','SPM-Print','HandleVisibility','off','CallBack','spm_figure(''Print'',gcbf)');
t1=uimenu(t0,'Label','default print file','HandleVisibility','off','CallBack','spm_figure(''Print'',gcbf)');
t1=uimenu(t0,'Label','other print file','HandleVisibility','off','CallBack','spm_figure(''PrintTo'',spm_figure(''FindWin'',''Graphics''))');


% ### CODE FOR SATELLITE FIGURE ###
% Code checks if there is a satellite window and if results are currently displayed
% It assumes that if hReg is invalid then there are no results currently displayed
% Modified by DRG to display a satellite figure 02/14/01.
cb = ['global SatWindow,',...
        'try,',...
        'tmp = get(hReg);,',...
        'ResFlag = 1;',...
        'catch,',...
        'ResFlag = 0;',...
        'end,',...
        'if SatWindow,',...
        'figure(SatWindow),',...
        'else,',...
        'if ResFlag,',...
        'spm_setup_satfig,',...
        'end,',...
        'end'];
uimenu( F,'Label','Results-Fig','HandleVisibility','off','Callback',cb);
% ### END NEW CODE ###


set(0,'ShowHiddenHandles',cSHH)
try
    spm_jobman('pulldown');
end

%=======================================================================


case 'figcontextmenu'
%=======================================================================
% h = spm_figure('FigContextMenu',F)
if nargin<2
    F = get(0,'CurrentFigure');
    if isempty(F), error('no figure'), end
else
    F = spm_figure('FindWin',varargin{2});
    if isempty(F), error('no such figure'), end
end
h       = uicontextmenu('Parent',F,'HandleVisibility','CallBack');
cSHH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')
copy_menu(F,h);
set(0,'ShowHiddenHandles',cSHH)
set(F,'UIContextMenu',h)
varargout = {h};


case 'colormap'
%=======================================================================
% spm_figure('ColorMap',ColAction,h)
if nargin<3, h=[]; else h=varargin{3}; end
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

case 'graphicshandle'
%=======================================================================
% h = spm_figure('GraphicsHandle',F)
if nargin<2, F=gcbf; else F=spm_figure('FindWin',varargin{2}); end
if isempty(F), return, end

tmp = get(F,'Name');
set(F,'Name',...
    'Handle: Select item to identify, MiddleMouse=parent, RightMouse=cancel...');
set(F,'Pointer','CrossHair')
waitforbuttonpress;
h        = gco(F);
hType    = get(h,'Type');
SelnType = get(gcf,'SelectionType');
set(F,'Pointer','Arrow','Name',tmp)

if ~strcmp(SelnType,'alt') && ~isempty(h) && gcf==F
    str = sprintf('Selected (%s) object',get(h,'Type'));
    if strcmp(SelnType,'normal')
        str = sprintf('%s: handle',str);
    else
        h = get(h,'Parent');
        str = sprintf('%s: handle of parent (%s) object',str,get(h,'Type'));
    end
    if nargout==0
        assignin('base','ans',h)
        fprintf('\n%s: \n',str)
        ans = h
    else
        varargout={h};
    end
else
    varargout={[]};
end



otherwise
%=======================================================================
warning(['Illegal Action string: ',Action])
end
return;
%=======================================================================


%=======================================================================
function copy_menu(F,G)
%=======================================================================
handles = findobj(get(F,'Children'),'Flat','Type','uimenu','Visible','on');
if isempty(handles), return; end;
for F1=handles',
    if ~strcmp(get(F1,'Label'),'&Window'),
        G1 = uimenu(G,'Label',get(F1,'Label'),...
            'CallBack',get(F1,'CallBack'),...
            'Position',get(F1,'Position'),...
            'Separator',get(F1,'Separator'));
        copy_menu(F1,G1);
    end;
end;
return;
%=======================================================================
