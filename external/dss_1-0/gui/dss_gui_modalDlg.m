function varargout = dss_gui_modalDlg(varargin)
% DSS_GUI_MODALDLG M-file for dss_gui_modalDlg.fig
%      DSS_GUI_MODALDLG by itself, creates a new DSS_GUI_MODALDLG or raises the
%      existing singleton*.
%
%      H = DSS_GUI_MODALDLG returns the handle to a new DSS_GUI_MODALDLG or the handle to
%      the existing singleton*.
%
%      DSS_GUI_MODALDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSS_GUI_MODALDLG.M with the given input arguments.
%
%      DSS_GUI_MODALDLG('Property','Value',...) creates a new DSS_GUI_MODALDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dss_gui_modalDlg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dss_gui_modalDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Includes functionality for a multi purpose Yes/No question window.
%      This file is used by other DSS files and should not be used alone.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% Edit the above text to modify the response to help dss_gui_modalDlg

% Last Modified by GUIDE v2.5 31-Jan-2005 10:19:09

% $Id: dss_gui_modalDlg.m,v 1.7 2005/04/20 10:19:24 kosti Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dss_gui_modalDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @dss_gui_modalDlg_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before dss_gui_modalDlg is made visible.
function dss_gui_modalDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dss_gui_modalDlg (see VARARGIN)

% Choose default command line output for dss_gui_modalDlg
handles.output = 'Yes';

% Update handles structure
guidata(hObject, handles);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
         case 'title'
          set(hObject, 'Name', varargin{index+1});
         case 'string'
          set(handles.text1, 'String', varargin{index+1});
        end
    end
end

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Show a question icon from dialogicons.mat - variables questIconData
% and questIconMap
load dialogicons.mat

IconData=questIconData;
questIconMap(256,:) = get(handles.modal_dlg, 'Color');
IconCMap=questIconMap;

Img=image(IconData, 'Parent', handles.axes1);
set(handles.modal_dlg, 'Colormap', IconCMap);

set(handles.axes1, ...
    'Visible', 'off', ...
    'YDir'   , 'reverse'       , ...
    'XLim'   , get(Img,'XData'), ...
    'YLim'   , get(Img,'YData')  ...
    );

% Make the GUI modal
set(handles.modal_dlg,'WindowStyle','modal')

% UIWAIT makes dss_gui_modalDlg wait for user response (see UIRESUME)
uiwait(handles.modal_dlg);

% --- Outputs from this function are returned to the command line.
function varargout = dss_gui_modalDlg_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.modal_dlg);


% --- Executes when user attempts to close figure1.
function modal_dlg_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.modal_dlg, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.modal_dlg);
else
    % The GUI is no longer waiting, just close it
    delete(handles.modal_dlg);
end


% --- Executes on key press over figure1 with no controls selected.
function modal_dlg_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.modal_dlg);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.modal_dlg);
end    


% --- Executes on button press in yesButton.
function yesButton_Callback(hObject, eventdata, handles)
% hObject    handle to yesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.modal_dlg);

% --- Executes on button press in noButton.
function noButton_Callback(hObject, eventdata, handles)
% hObject    handle to noButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.modal_dlg);