function varargout = dss_gui_insertFunc(varargin)
% DSS_GUI_INSERTFUNC M-file for dss_gui_insertFunc.fig
%      DSS_GUI_INSERTFUNC, by itself, creates a new DSS_GUI_INSERTFUNC or raises the existing
%      singleton*.
%
%      H = DSS_GUI_INSERTFUNC returns the handle to a new DSS_GUI_INSERTFUNC or the handle to
%      the existing singleton*.
%
%      DSS_GUI_INSERTFUNC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSS_GUI_INSERTFUNC.M with the given input arguments.
%
%      DSS_GUI_INSERTFUNC('Property','Value',...) creates a new DSS_GUI_INSERTFUNC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dss_gui_insertFunc_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dss_gui_insertFunc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Includes functionality of the insert function window.
%      This file is used by other DSS files and should not be used alone.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% Edit the above text to modify the response to help dss_gui_insertFunc

% Last Modified by GUIDE v2.5 25-Jan-2005 14:42:03

% $Id: dss_gui_insertFunc.m,v 1.8 2005/04/20 10:19:24 kosti Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dss_gui_insertFunc_OpeningFcn, ...
                   'gui_OutputFcn',  @dss_gui_insertFunc_OutputFcn, ...
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


% --- Executes just before dss_gui_insertFunc is made visible.
function dss_gui_insertFunc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dss_gui_insertFunc (see VARARGIN)

% Choose default command line output for dss_gui_insertFunc
handles.output = hObject;

% Name dialog
set(hObject, 'Name', 'Custom function');

% Handles from main GUI
handles.mainGUIhandles = varargin{2};
handles.caller = varargin{1};

% Find out who called this dialog
guidata(hObject, handles);
handles.callerTag = get(handles.caller, 'Tag');

% Make dialog modal
set(hObject, 'WindowStyle','modal');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dss_gui_insertFunc wait for user response (see UIRESUME)
% uiwait(handles.insertfGUI);


% --- Outputs from this function are returned to the command line.
function varargout = dss_gui_insertFunc_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function functionNameEditText_Callback(hObject, eventdata, handles)
% hObject    handle to functionNameEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of functionNameEditText as text
%        str2double(get(hObject,'String')) returns contents of functionNameEditText as a double


% --- Executes during object creation, after setting all properties.
function functionNameEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to functionNameEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DENOISING_FUNCTIONS_SORTED;
global ORTHO_FUNCTIONS_SORTED;
global PREPRO_FUNCTIONS_SORTED;
global ALPHA_FUNCTIONS_SORTED;
global BETA_FUNCTIONS_SORTED;
global GAMMA_FUNCTIONS_SORTED;
global CUSTOM_ALPHA_FUNCTIONS;
global CUSTOM_BETA_FUNCTIONS;
global CUSTOM_GAMMA_FUNCTIONS;

% Get users input
func_name = get(handles.functionNameEditText, 'String');
% Make sure it's not empty input
if ~isempty(func_name)
    % Check if file exists
   if exist(func_name) == 2 
      % Remove file extension (.m) if inserted
      end_char = length(func_name);
      if strcmp(func_name(end_char), 'm') & strcmp(func_name(end_char-1), '.')
          func_name = func_name(1:end_char-2);
      end
      % Who is the caller
      if strcmp(get(handles.caller, 'Tag'), 'selectDenoisingPopup')
            % selectDenoisingPopup called
            strings = get(handles.mainGUIhandles.selectDenoisingPopup, 'String');
            numberOfStrings = size(strings, 1);
            strings{numberOfStrings, 1} = func_name;
            strings{numberOfStrings+1, 1} = 'Custom';
            set(handles.mainGUIhandles.selectDenoisingPopup, 'String', strings);
            set(handles.mainGUIhandles.selectDenoisingPopup, 'Value', numberOfStrings);
            DENOISING_FUNCTIONS_SORTED{length(DENOISING_FUNCTIONS_SORTED)+1} = func_name;
            if isfield(handles.mainGUIhandles, 'customDenoisingFunctions')
                n = length(handles.mainGUIhandles.customDenoisingFunctions);
                handles.mainGUIhandles.customDenoisingFunctions{n+1} = func_name;
            else
                handles.mainGUIhandles.customDenoisingFunctions = {};
                handles.mainGUIhandles.customDenoisingFunctions{1} = func_name;
            end
            guidata(handles.mainGUIhandles.mainGUI, handles.mainGUIhandles);
      elseif strcmp(get(handles.caller, 'Tag'), 'selectAlphafPopup')
            % selectAlphafPopup called
            strings = get(handles.mainGUIhandles.selectAlphafPopup, 'String');
            numberOfStrings = size(strings, 1);
            strings{numberOfStrings, 1} = func_name;
            strings{numberOfStrings+1, 1} = 'Custom';
            set(handles.mainGUIhandles.selectAlphafPopup, 'String', strings);
            set(handles.mainGUIhandles.selectAlphafPopup, 'Value', numberOfStrings);
            ALPHA_FUNCTIONS_SORTED{length(ALPHA_FUNCTIONS_SORTED)+1} = func_name;
            if isempty(CUSTOM_ALPHA_FUNCTIONS)
                CUSTOM_ALPHA_FUNCTIONS = {};
                CUSTOM_ALPHA_FUNCTIONS{1} = func_name;
            else
                CUSTOM_ALPHA_FUNCTIONS{length(CUSTOM_ALPHA_FUNCTIONS)+1} = func_name;
            end
      elseif strcmp(get(handles.caller, 'Tag'), 'selectBetafPopup')
            % selectBetafPopup called
            strings = get(handles.mainGUIhandles.selectBetafPopup, 'String');
            numberOfStrings = size(strings, 1);
            strings{numberOfStrings, 1} = func_name;
            strings{numberOfStrings+1, 1} = 'Custom';
            set(handles.mainGUIhandles.selectBetafPopup, 'String', strings);
            set(handles.mainGUIhandles.selectBetafPopup, 'Value', numberOfStrings);
            BETA_FUNCTIONS_SORTED{length(BETA_FUNCTIONS_SORTED)+1} = func_name;
            if isempty(CUSTOM_BETA_FUNCTIONS)
                CUSTOM_BETA_FUNCTIONS = {};
                CUSTOM_BETA_FUNCTIONS{1} = func_name;
            else
                CUSTOM_BETA_FUNCTIONS{length(CUSTOM_BETA_FUNCTIONS)+1} = func_name;
            end
      elseif strcmp(get(handles.caller, 'Tag'), 'selectGammafPopup')
            % selectGammafPopup called
            strings = get(handles.mainGUIhandles.selectGammafPopup, 'String');
            numberOfStrings = size(strings, 1);
            strings{numberOfStrings, 1} = func_name;
            strings{numberOfStrings+1, 1} = 'Custom';
            set(handles.mainGUIhandles.selectGammafPopup, 'String', strings);
            set(handles.mainGUIhandles.selectGammafPopup, 'Value', numberOfStrings);
            GAMMA_FUNCTIONS_SORTED{length(GAMMA_FUNCTIONS_SORTED)+1} = func_name;
            if isempty(CUSTOM_GAMMA_FUNCTIONS)
                CUSTOM_GAMMA_FUNCTIONS = {};
                CUSTOM_GAMMA_FUNCTIONS{1} = func_name;
            else
                CUSTOM_GAMMA_FUNCTIONS{length(CUSTOM_GAMMA_FUNCTIONS)+1} = func_name;
            end
      elseif strcmp(get(handles.caller, 'Tag'), 'selectOrthoPopup')
            % selectOrthoPopup called
            strings = get(handles.mainGUIhandles.selectOrthoPopup, 'String');
            numberOfStrings = size(strings, 1);
            strings{numberOfStrings, 1} = func_name;
            strings{numberOfStrings+1, 1} = 'Custom';
            set(handles.mainGUIhandles.selectOrthoPopup, 'String', strings);
            set(handles.mainGUIhandles.selectOrthoPopup, 'Value', numberOfStrings);  
            ORTHO_FUNCTIONS_SORTED{length(ORTHO_FUNCTIONS_SORTED)+1} = func_name;
      elseif strcmp(get(handles.caller, 'Tag'), 'selectWhiteningPopup')
            % selectWhiteningPopup called
            strings = get(handles.mainGUIhandles.selectWhiteningPopup, 'String');
            numberOfStrings = size(strings, 1);
            strings{numberOfStrings, 1} = func_name;
            strings{numberOfStrings+1, 1} = 'Custom';
            set(handles.mainGUIhandles.selectWhiteningPopup, 'String', strings);
            set(handles.mainGUIhandles.selectWhiteningPopup, 'Value', numberOfStrings);
            PREPRO_FUNCTIONS_SORTED{length(PREPRO_FUNCTIONS_SORTED)+1} = func_name;
      end
      % Close dialog
      close(handles.insertfGUI);
   else
       errordlg('File does not exist.', 'Input Error', 'modal');
   end
else
    errordlg('Empty input.', 'Input Error', 'modal');
end

% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LAST_DENOISING_FUNC;
global LAST_WHITENING_FUNC;
global LAST_ORTHO_FUNC;

% Set default item as selected item, so 'Custom' won't stay selected
if strcmp(get(handles.caller, 'Tag'), 'selectWhiteningPopup')
    set(handles.mainGUIhandles.selectWhiteningPopup, 'Value', LAST_WHITENING_FUNC);
elseif strcmp(get(handles.caller, 'Tag'), 'selectDenoisingPopup')
    set(handles.mainGUIhandles.selectDenoisingPopup, 'Value', LAST_DENOISING_FUNC);
elseif strcmp(get(handles.caller, 'Tag'), 'selectOrthoPopup')
    set(handles.mainGUIhandles.selectOrthoPopup, 'Value', LAST_ORTHO_FUNC);
elseif strcmp(get(handles.caller, 'Tag'), 'selectAlphafPopup')
    set(handles.mainGUIhandles.selectAlphafPopup, 'Value', 1);
    set(handles.mainGUIhandles.fixedAlphaEditText, 'Enable', 'on', 'Visible', 'on');
elseif strcmp(get(handles.caller, 'Tag'), 'selectBetafPopup')
    set(handles.mainGUIhandles.selectBetafPopup, 'Value', 1);
    set(handles.mainGUIhandles.fixedBetaEditText, 'Enable', 'on', 'Visible', 'on');
elseif strcmp(get(handles.caller, 'Tag'), 'selectGammafPopup')
    set(handles.mainGUIhandles.selectGammafPopup, 'Value', 1);
    set(handles.mainGUIhandles.fixedGammaEditText, 'Enable', 'on', 'Visible', 'on');
end
    
% Close dialog
close(handles.insertfGUI);