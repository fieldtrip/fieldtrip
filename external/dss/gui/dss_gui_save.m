function varargout = dss_gui_save(varargin)
% DSS_GUI_SAVE M-file for dss_gui_save.fig
%      DSS_GUI_SAVE, by itself, creates a new DSS_GUI_SAVE or raises the existing
%      singleton*.
%
%      H = DSS_GUI_SAVE returns the handle to a new DSS_GUI_SAVE or the handle to
%      the existing singleton*.
%
%      DSS_GUI_SAVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSS_GUI_SAVE.M with the given input arguments.
%
%      DSS_GUI_SAVE('Property','Value',...) creates a new DSS_GUI_SAVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dss_gui_save_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dss_gui_save_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Includes functionality for the save results window.
%      This file is used by other DSS files and should not be used alone.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% Edit the above text to modify the response to help dss_gui_save

% Last Modified by GUIDE v2.5 03-Mar-2005 09:25:21

% $Id: dss_gui_save.m,v 1.8 2005/04/20 10:19:24 kosti Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dss_gui_save_OpeningFcn, ...
                   'gui_OutputFcn',  @dss_gui_save_OutputFcn, ...
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


% --- Executes just before dss_gui_save is made visible.
function dss_gui_save_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dss_gui_save (see VARARGIN)

% Choose default command line output for dss_gui_save
handles.output = hObject;

% Name dialog
set(hObject, 'Name', 'Save results');

% Handles from main GUI
handles.mainGUIhandles = varargin{1};

% Make dialog modal
set(hObject, 'WindowStyle','modal');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dss_gui_save wait for user response (see UIRESUME)
% uiwait(handles.saveGUI);


% --- Outputs from this function are returned to the command line.
function varargout = dss_gui_save_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function suffixEditText_Callback(hObject, eventdata, handles)
% hObject    handle to suffixEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of suffixEditText as text
%        str2double(get(hObject,'String')) returns contents of suffixEditText as a double


% --- Executes during object creation, after setting all properties.
function suffixEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to suffixEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inputString = get(handles.suffixEditText, 'String');
string_A = strcat('A_', inputString);
string_W = strcat('W_', inputString);
string_S = strcat('S_', inputString);
string_state = strcat('state_', inputString);
A_duplicate = false;
W_duplicate = false;
S_duplicate = false;
state_duplicate = false;

% Use default names in case of an empty input
if isempty(inputString)
    string_A = 'A_DSS';
    string_W = 'W_DSS';
    string_S = 'S_DSS';
    string_state = 'state_DSS';
end

    % Check if there are already variables by the same name
    vars = evalin('base', 'who');
    vars_length = length(vars);
    for i = 1:vars_length
        if strcmp(string_A, vars(i))
            A_duplicate = true;
        end
        if strcmp(string_W, vars(i))
            W_duplicate = true;
        end
        if strcmp(string_S, vars(i))
            S_duplicate = true;
        end
        if strcmp(string_state, vars(i))
            state_duplicate = true;
        end
    end
    
    if A_duplicate
        % Ask for overwrite
        response = dss_gui_modalDlg('Title', ['Confirm Overwrite: ' string_A]);
        switch response
            case 'Yes'
                var = handles.mainGUIhandles.A;
                assignin('base', string_A, var);
            case 'No'
                % Do nothing
        end
    else
        % Write variable
        var = handles.mainGUIhandles.A;
        assignin('base', string_A, var);
    end
    if W_duplicate
        % Ask for overwrite
        response = dss_gui_modalDlg('Title', ['Confirm Overwrite: ' string_W]);
        switch response
            case 'Yes'
                var = handles.mainGUIhandles.W;
                assignin('base', string_W, var);
            case 'No'
                % Do nothing
        end
    else
        % Write variable
        var = handles.mainGUIhandles.W;
        assignin('base', string_W, var);
    end
    if S_duplicate
        % Ask for overwrite
        response = dss_gui_modalDlg('Title', ['Confirm Overwrite: ' string_S]);
        switch response
            case 'Yes'
                var = handles.mainGUIhandles.S;
                assignin('base', string_S, var);
            case 'No'
                % Do nothing
        end
    else
        % Write variable
        var = handles.mainGUIhandles.S;
        assignin('base', string_S, var);
    end
    if state_duplicate
        % Ask for overwrite
        response = dss_gui_modalDlg('Title', ['Confirm Overwrite: ' string_state]);
        switch response
            case 'Yes'
                var = handles.mainGUIhandles.state;
                assignin('base', string_state, var);
            case 'No'
                % Do nothing
        end
    else
        % Write variable
        var = handles.mainGUIhandles.state;
        assignin('base', string_state, var);
    end

% Close dialog
close(handles.saveGUI);

% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.saveGUI);


% --- Executes on button press in plotButton.
function plotButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reporting
fig = figure; clf;
set(fig, 'Name', 'DSS Results');
try
    test_report_result(handles.mainGUIhandles.state.S, handles.mainGUIhandles.state);
catch
    errordlg(lasterr, 'DSS Report Function Error', 'modal');
end