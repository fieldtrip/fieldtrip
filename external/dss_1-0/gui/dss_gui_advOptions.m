function varargout = dss_gui_advOptions(varargin)
% DSS_GUI_ADVOPTIONS M-file for dss_gui_advOptions.fig
%      DSS_GUI_ADVOPTIONS, by itself, creates a new DSS_GUI_ADVOPTIONS or raises the existing
%      singleton*.
%
%      H = DSS_GUI_ADVOPTIONS returns the handle to a new DSS_GUI_ADVOPTIONS or the handle to
%      the existing singleton*.
%
%      DSS_GUI_ADVOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSS_GUI_ADVOPTIONS.M with the given input arguments.
%
%      DSS_GUI_ADVOPTIONS('Property','Value',...) creates a new DSS_GUI_ADVOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dss_gui_advOptions_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dss_gui_advOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Includes functionality of the advanced options window.
%      This file is used by other DSS files and should not be used alone.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% Edit the above text to modify the response to help dss_gui_advOptions

% Last Modified by GUIDE v2.5 28-Jan-2005 15:54:48

% $Id: dss_gui_advOptions.m,v 1.8 2005/04/20 10:19:24 kosti Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dss_gui_advOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @dss_gui_advOptions_OutputFcn, ...
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


% --- Executes just before dss_gui_advOptions is made visible.
function dss_gui_advOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dss_gui_advOptions (see VARARGIN)

% Choose default command line output for dss_gui_advOptions
handles.output = hObject;

global ALPHA_FUNCTIONS_SORTED;
global BETA_FUNCTIONS_SORTED;
global GAMMA_FUNCTIONS_SORTED;

% Set dialog name
set(hObject, 'Name', 'Advanced options');

% Make dialog modal
set(hObject, 'WindowStyle','modal');

% Handles from mainGUI
handles.mainGUI = varargin{1};
handles.mainGUIhandles = varargin{2};

% Populate popUp menus by reading all available m-files
cleanPopup(handles.selectAlphafPopup);
cleanPopup(handles.selectBetafPopup);
cleanPopup(handles.selectGammafPopup);
readFunctions(handles, 'alpha');
readFunctions(handles, 'beta');
readFunctions(handles, 'gamma');

% Flag for monitoring modifications
handles.modified_flag = false;

% Update handles structure
guidata(hObject, handles);

% Reading current parameters

% Alpha function
if isfield(handles.mainGUIhandles.params, 'alphaf')
    func_name = func2str(handles.mainGUIhandles.params.alphaf.h);
    % Check if the function is already in the list
    found = 0;
    for i = 1:length(ALPHA_FUNCTIONS_SORTED)
        if strcmp(func_name, ALPHA_FUNCTIONS_SORTED{i})
            set(handles.selectAlphafPopup, 'Value', i+1);
            found = 1;
        end
    end
    if ~found
        % Update custom function to the list
        cleanPopup(handles.selectAlphafPopup);
        itemToPopup(handles.selectAlphafPopup, func_name);
    end
    
    % Disable fixed alpha
    set(handles.fixedAlphaEditText, 'Enable', 'off', 'Visible', 'off');
else
    % Enable fixed alpha
    set(handles.fixedAlphaEditText, 'Enable', 'on', 'Visible', 'on');
    set(handles.selectAlphafPopup, 'Value', 1);
end

% Beta function
if isfield(handles.mainGUIhandles.params, 'betaf')
    func_name = func2str(handles.mainGUIhandles.params.betaf.h);
    % Check if the function is already in the list
    found = 0;
    for i = 1:length(BETA_FUNCTIONS_SORTED)
        if strcmp(func_name, BETA_FUNCTIONS_SORTED{i})
            set(handles.selectBetafPopup, 'Value', i+1);
            found = 1;
        end
    end
    if ~found
        % Update custom function to the list
        cleanPopup(handles.selectBetafPopup);
        itemToPopup(handles.selectBetafPopup, func_name);
    end

    % Disable fixed beta
    set(handles.fixedBetaEditText, 'Enable', 'off', 'Visible', 'off');
else
    % Enable fixed beta
    set(handles.fixedBetaEditText, 'Enable', 'on', 'Visible', 'on');
    set(handles.selectBetafPopup, 'Value', 1);
end

% Gamma function
if isfield(handles.mainGUIhandles.params, 'gammaf')
    func_name = func2str(handles.mainGUIhandles.params.gammaf.h);
    % Check if the function is already in the list
    found = 0;
    for i = 1:length(GAMMA_FUNCTIONS_SORTED)
        if strcmp(func_name, GAMMA_FUNCTIONS_SORTED{i})
            set(handles.selectGammafPopup, 'Value', i+1);
            found = 1;
        end
    end
    if ~found
        % Update custom function to the list
        cleanPopup(handles.selectGammafPopup);
        itemToPopup(handles.selectGammafPopup, func_name);
    end

    % Disable fixed gamma
    set(handles.fixedGammaEditText, 'Enable', 'off', 'Visible', 'off');
else
    % Enable fixed gamma
    set(handles.fixedGammaEditText, 'Enable', 'on', 'Visible', 'on');
    set(handles.selectGammafPopup, 'Value', 1);
end

% Fixed alpha, beta, gamma
if isfield(handles.mainGUIhandles.params, 'alpha')
    set(handles.fixedAlphaEditText, 'String', num2str(handles.mainGUIhandles.params.alpha));
end
if isfield(handles.mainGUIhandles.params, 'beta')
    set(handles.fixedBetaEditText, 'String', num2str(handles.mainGUIhandles.params.beta));
end
if isfield(handles.mainGUIhandles.params, 'gamma')
    set(handles.fixedGammaEditText, 'String', num2str(handles.mainGUIhandles.params.gamma));
end

% Verbosity
if isfield(handles.mainGUIhandles.params, 'verbose')
    set(handles.verbosityPopup, 'Value', handles.mainGUIhandles.params.verbose + 1);
else
    set(handles.verbosityPopup, 'Value', 1);
end

% Interruptability
if isfield(handles.mainGUIhandles.params, 'gui_interrupt')
    set(handles.interruptableToggle, 'Value', handles.mainGUIhandles.params.gui_interrupt);
end

% UIWAIT makes dss_gui_advOptions wait for user response (see UIRESUME)
% uiwait(handles.advOptGUI);


% --- Outputs from this function are returned to the command line.
function varargout = dss_gui_advOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in selectAlphafPopup.
function selectAlphafPopup_Callback(hObject, eventdata, handles)
% hObject    handle to selectAlphafPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectAlphafPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectAlphafPopup

% Enable apply button
set(handles.applyButton, 'Enable', 'on');
handles.modified_flag = true;
guidata(hObject, handles);

% Call insertFunction if 'Custom' is selected
alphaItems = get(hObject, 'String');
index_selected = get(hObject, 'Value');
if strcmp(alphaItems{index_selected(1)}, 'Custom')
    dss_gui_insertFunc(hObject, handles);
end

% Enable and disable optional elements
if (get(hObject, 'Value')) ~= 1
    % Disable fixed alpha
    set(handles.fixedAlphaEditText, 'Enable', 'off', 'Visible', 'off');
else
    % Enable fixed alpha
    set(handles.fixedAlphaEditText, 'Enable', 'on', 'Visible', 'on');
end

% --- Executes during object creation, after setting all properties.
function selectAlphafPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectAlphafPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in selectBetafPopup.
function selectBetafPopup_Callback(hObject, eventdata, handles)
% hObject    handle to selectBetafPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectBetafPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectBetafPopup

% Enable apply button
set(handles.applyButton, 'Enable', 'on');
handles.modified_flag = true;
guidata(hObject, handles);

% Call insertFunc if 'Custom' is selected
betaItems = get(hObject, 'String');
index_selected = get(hObject, 'Value');
if strcmp(betaItems{index_selected(1)}, 'Custom')
    dss_gui_insertFunc(hObject, handles);
end

% Enable and disable optional elements
if (get(hObject, 'Value')) ~= 1
    % Disable fixed beta
    set(handles.fixedBetaEditText, 'Enable', 'off', 'Visible', 'off');
else
    % Enable fixed beta
    set(handles.fixedBetaEditText, 'Enable', 'on', 'Visible', 'on');
end

% --- Executes during object creation, after setting all properties.
function selectBetafPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectBetafPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in selectGammafPopup.
function selectGammafPopup_Callback(hObject, eventdata, handles)
% hObject    handle to selectGammafPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectGammafPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectGammafPopup

% Enable apply button
set(handles.applyButton, 'Enable', 'on');
handles.modified_flag = true;
guidata(hObject, handles);

% Call insertFunction if 'Custom' is selected
gammaItems = get(hObject, 'String');
index_selected = get(hObject, 'Value');
if strcmp(gammaItems{index_selected(1)}, 'Custom')
    dss_gui_insertFunc(hObject, handles);
end

% Enable and disable optional elements
if (get(hObject, 'Value')) ~= 1
    % Disable fixed gamma
    set(handles.fixedGammaEditText, 'Enable', 'off', 'Visible', 'off');
else
    % Enable fixed gamma
    set(handles.fixedGammaEditText, 'Enable', 'on', 'Visible', 'on');
end


% --- Executes during object creation, after setting all properties.
function selectGammafPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectGammafPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in interruptableToggle.
function interruptableToggle_Callback(hObject, eventdata, handles)
% hObject    handle to interruptableToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of interruptableToggle
set(handles.applyButton, 'Enable', 'on');
handles.modified_flag = true;
guidata(hObject, handles);


function fixedAlphaEditText_Callback(hObject, eventdata, handles)
% hObject    handle to fixedAlphaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixedAlphaEditText as text
%        str2double(get(hObject,'String')) returns contents of fixedAlphaEditText as a double

% Enable apply button
set(handles.applyButton, 'Enable', 'on');
handles.modified_flag = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fixedAlphaEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedAlphaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function fixedBetaEditText_Callback(hObject, eventdata, handles)
% hObject    handle to fixedBetaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixedBetaEditText as text
%        str2double(get(hObject,'String')) returns contents of fixedBetaEditText as a double

% Enable apply button
set(handles.applyButton, 'Enable', 'on');
handles.modified_flag = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fixedBetaEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedBetaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function fixedGammaEditText_Callback(hObject, eventdata, handles)
% hObject    handle to fixedGammaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixedGammaEditText as text
%        str2double(get(hObject,'String')) returns contents of fixedGammaEditText as a double

% Enable apply button
set(handles.applyButton, 'Enable', 'on');
handles.modified_flag = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fixedGammaEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedGammaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Something has been modified, inform mainGUI
if handles.modified_flag
    handles.mainGUIhandles.paramDefaultsFlag = false;
    set(handles.mainGUIhandles.paramInfoText, 'String', 'Custom');
    set(handles.mainGUIhandles.paramResetButton, 'Enable', 'on');
    guidata(handles.mainGUIhandles.mainGUI, handles.mainGUIhandles);
    handles.modified_flag = false;
    guidata(hObject, handles);
end

loadData(handles);

% Close dialog
close(handles.advOptGUI);


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.advOptGUI);


% --- Executes on button press in defaultsButton.
function defaultsButton_Callback(hObject, eventdata, handles)
% hObject    handle to defaultsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALPHA_FUNCTIONS_SORTED;
global BETA_FUNCTIONS_SORTED;
global GAMMA_FUNCTIONS_SORTED;
global CUSTOM_ALPHA_FUNCTIONS;
global CUSTOM_BETA_FUNCTIONS;
global CUSTOM_GAMMA_FUNCTIONS;

% Populate popUp menus by reading all available m-files
cleanPopup(handles.selectAlphafPopup);
cleanPopup(handles.selectBetafPopup);
cleanPopup(handles.selectGammafPopup);
readFunctions(handles, 'alpha');
readFunctions(handles, 'beta');
readFunctions(handles, 'gamma');

% Disable apply button
set(handles.applyButton, 'Enable', 'off');
handles.modified_flag = false;
guidata(hObject, handles);

if handles.mainGUIhandles.state_given
    % Set default values from given state
    % ALPHA, BETA, GAMMA
    % Functions
    % Alphaf
    found = false;
    if isfield(handles.mainGUIhandles.state, 'alphaf')
        if ~isempty(handles.mainGUIhandles.state.alphaf)
            handles.mainGUIhandles.params.alphaf = handles.mainGUIhandles.state.alphaf;
            % Select proper function from popUp
            if isfield(handles.mainGUIhandles.state.alphaf, 'h')
                for i = 1:length(ALPHA_FUNCTIONS_SORTED)
                    if strcmp(ALPHA_FUNCTIONS_SORTED{i}, func2str(handles.mainGUIhandles.state.alphaf.h))
                        set(handles.selectAlphafPopup, 'Value', i+1);
                        found = true;
                        break;
                    end
                end
                if ~found
                    ALPHA_FUNCTIONS_SORTED{length(ALPHA_FUNCTIONS_SORTED)+1} = func2str(handles.mainGUIhandles.state.alphaf.h);
                    itemToPopup(handles.selectAlphafPopup, func2str(handles.mainGUIhandles.state.alphaf.h));
                end
            end
            % Disable fixed alpha
            set(handles.fixedAlphaEditText, 'Enable', 'off', 'Visible', 'off');
        else
            if isfield(handles.mainGUIhandles.params, 'alphaf')
                handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'alphaf');
            end
            % Enable fixed alpha
            set(handles.fixedAlphaEditText, 'Enable', 'on', 'Visible', 'on');
            set(handles.selectAlphafPopup, 'Value', 1);
        end
    else
        if isfield(handles.mainGUIhandles.params, 'alphaf')
            handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'alphaf');
        end
        % Enable fixed alpha
        set(handles.fixedAlphaEditText, 'Enable', 'on', 'Visible', 'on');
        set(handles.selectAlphafPopup, 'Value', 1);
    end
    found = false;
    % Betaf
    if isfield(handles.mainGUIhandles.state, 'betaf')
        if ~isempty(handles.mainGUIhandles.state.betaf)
            handles.mainGUIhandles.params.betaf = handles.mainGUIhandles.state.betaf;
            % Select proper function from popUp
            if isfield(handles.mainGUIhandles.state.betaf, 'h')
                for i = 1:length(BETA_FUNCTIONS_SORTED)
                    if strcmp(BETA_FUNCTIONS_SORTED{i}, func2str(handles.mainGUIhandles.state.betaf.h))
                        set(handles.selectBetafPopup, 'Value', i+1);
                        found = true;
                        break;
                    end
                end
                if ~found
                    BETA_FUNCTIONS_SORTED{length(BETA_FUNCTIONS_SORTED)+1} = func2str(handles.mainGUIhandles.state.betaf.h);
                    itemToPopup(handles.selectBetafPopup, func2str(handles.mainGUIhandles.state.betaf.h));
                end
            end
            % Disable fixed beta
            set(handles.fixedBetaEditText, 'Enable', 'off', 'Visible', 'off');
        else
            if isfield(handles.mainGUIhandles.params, 'betaf')
                handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'betaf');
            end
            % Enable fixed beta
            set(handles.fixedBetaEditText, 'Enable', 'on', 'Visible', 'on');
            set(handles.selectBetafPopup, 'Value', 1);
        end
    else
        if isfield(handles.mainGUIhandles.params, 'betaf')
            handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'betaf');
        end
        % Enable fixed beta
        set(handles.fixedBetaEditText, 'Enable', 'on', 'Visible', 'on');
        set(handles.selectBetafPopup, 'Value', 1);
    end
    found = false;
    % Gammaf
    if isfield(handles.mainGUIhandles.state, 'gammaf')
        if ~isempty(handles.mainGUIhandles.state.gammaf)
            handles.mainGUIhandles.params.gammaf = handles.mainGUIhandles.state.gammaf;
            % Select proper function from popUp
            functions = get(handles.selectGammafPopup, 'String');
            if isfield(handles.mainGUIhandles.state.gammaf, 'h')
                for i = 1:length(GAMMA_FUNCTIONS_SORTED)
                    if strcmp(GAMMA_FUNCTIONS_SORTED{i}, func2str(handles.mainGUIhandles.state.gammaf.h))
                        set(handles.selectGammafPopup, 'Value', i+1);
                        found = true;
                        break;
                    end
                end 
                if ~found
                    GAMMA_FUNCTIONS_SORTED{length(GAMMA_FUNCTIONS_SORTED)+1} = func2str(handles.mainGUIhandles.state.gammaf.h);
                    itemToPopup(handles.selectGammafPopup, func2str(handles.mainGUIhandles.state.gammaf.h));
                end
            end
            % Disable fixed gamma
            set(handles.fixedGammaEditText, 'Enable', 'off', 'Visible', 'off');
        else
            if isfield(handles.mainGUIhandles.params, 'gammaf')
                handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'gammaf');
            end
            % Enable fixed gamma
            set(handles.fixedGammaEditText, 'Enable', 'on', 'Visible', 'on');
            set(handles.selectGammafPopup, 'Value', 1);
        end
    else
        if isfield(handles.mainGUIhandles.params, 'gammaf')
            handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'gammaf');
        end
        % Enable fixed gamma
        set(handles.fixedGammaEditText, 'Enable', 'on', 'Visible', 'on');
        set(handles.selectGammafPopup, 'Value', 1);
    end
    
    % Fixed values
    if isfield(handles.mainGUIhandles.state, 'alpha')
        if ~isempty(handles.mainGUIhandles.state.alpha)
            handles.mainGUIhandles.params.alpha = handles.mainGUIhandles.state.alpha;
            set(handles.fixedAlphaEditText, 'String', num2str(handles.mainGUIhandles.state.alpha));
        end
    end
    if isfield(handles.mainGUIhandles.state, 'beta')
        if ~isempty(handles.mainGUIhandles.state.beta)
            handles.mainGUIhandles.params.beta = handles.mainGUIhandles.state.beta;
            set(handles.fixedBetaEditText, 'String', num2str(handles.mainGUIhandles.state.beta));
        end
    end
    if isfield(handles.mainGUIhandles.state, 'gamma')
        if ~isempty(handles.mainGUIhandles.state.gamma)
            handles.mainGUIhandles.params.gamma = handles.mainGUIhandles.state.gamma;
            set(handles.fixedGammaEditText, 'String', num2str(handles.mainGUIhandles.state.gamma));
        end
    end
    
    guidata(handles.mainGUIhandles.mainGUI, handles.mainGUIhandles);
else
    % Set default values
    set(handles.verbosityPopup, 'Value', 1);
    set(handles.interruptableToggle, 'Value', true);
    set(handles.selectAlphafPopup, 'Value', 1);
    set(handles.selectBetafPopup, 'Value', 1);
    set(handles.selectGammafPopup, 'Value', 1);
    set(handles.fixedAlphaEditText, 'String', '1');
    set(handles.fixedBetaEditText, 'String', '0');
    set(handles.fixedGammaEditText, 'String', '1');
    % Enable fixed alpha
    set(handles.fixedAlphaEditText, 'Enable', 'on', 'Visible', 'on');
    % Enable fixed beta
    set(handles.fixedBetaEditText, 'Enable', 'on', 'Visible', 'on');
    % Enable fixed gamma
    set(handles.fixedGammaEditText, 'Enable', 'on', 'Visible', 'on');
    % Start using the default values
    loadData(handles);
end


% --- Executes on selection change in verbosityPopup.
function verbosityPopup_Callback(hObject, eventdata, handles)
% hObject    handle to verbosityPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns verbosityPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from verbosityPopup
set(handles.applyButton, 'Enable', 'on');
handles.modified_flag = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function verbosityPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verbosityPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in applyButton.
function applyButton_Callback(hObject, eventdata, handles)
% hObject    handle to applyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject, 'Enable'), 'on')
    % Inform mainGUI if something has been modified
    if handles.modified_flag
        handles.mainGUIhandles.paramDefaultsFlag = false;
        set(handles.mainGUIhandles.paramInfoText, 'String', 'Custom');
        set(handles.mainGUIhandles.paramResetButton, 'Enable', 'on');
        guidata(handles.mainGUIhandles.mainGUI, handles.mainGUIhandles);
    end
    loadData(handles);
    set(hObject, 'Enable', 'off');
    handles.modified_flag = false;
    guidata(hObject, handles);
end


% --- Used with okButton and applyButton
function loadData(handles)

global ALPHA_FUNCTIONS_SORTED;
global BETA_FUNCTIONS_SORTED;
global GAMMA_FUNCTIONS_SORTED;

% Check for valid input and update mainGUI handles

% Alpha function
alphafunctions = get(handles.selectAlphafPopup, 'String');
index_selected = get(handles.selectAlphafPopup, 'Value');
switch alphafunctions{index_selected(1)}
    case 'Fixed value'
        if isfield(handles.mainGUIhandles.params, 'alphaf')
            handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'alphaf');
        end
    case 'Custom'
        errordlg('Invalid alpha function selection.', 'Input Error', 'modal');
        return;
    otherwise
        handles.mainGUIhandles.params.alphaf = struct('h', str2func(ALPHA_FUNCTIONS_SORTED{index_selected-1}));
end

% Beta function
betafunctions = get(handles.selectBetafPopup, 'String');
index_selected = get(handles.selectBetafPopup, 'Value');
switch betafunctions{index_selected(1)}
    case 'Fixed value'
        if isfield(handles.mainGUIhandles.params, 'betaf')
            handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'betaf');
        end
    case 'Custom'
        errordlg('Invalid beta function selection.', 'Input Error', 'modal');
        return;
    otherwise
        handles.mainGUIhandles.params.betaf = struct('h', str2func(BETA_FUNCTIONS_SORTED{index_selected-1}));
end

% Gamma function
gammafunctions = get(handles.selectGammafPopup, 'String');
index_selected = get(handles.selectGammafPopup, 'Value');
switch gammafunctions{index_selected(1)}
    case 'Fixed value'
        if isfield(handles.mainGUIhandles.params, 'gammaf')
            handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'gammaf');
        end
    case 'Custom'
        errordlg('Invalid gamma function selection.', 'Input Error', 'modal');
        return;
    otherwise
        handles.mainGUIhandles.params.gammaf = struct('h', str2func(GAMMA_FUNCTIONS_SORTED{index_selected-1}));
end

% Fixed alpha, beta , gamma
if strcmp(get(handles.fixedAlphaEditText, 'Enable'), 'on')
    handles.mainGUIhandles.params.alpha = str2num(get(handles.fixedAlphaEditText, 'String'));
else
    handles.mainGUIhandles.params.alpha = 1;
end
if strcmp(get(handles.fixedBetaEditText, 'Enable'), 'on')
    handles.mainGUIhandles.params.beta = str2num(get(handles.fixedBetaEditText, 'String'));
else
    handles.mainGUIhandles.params.beta = 0;
end
if strcmp(get(handles.fixedGammaEditText, 'Enable'), 'on')
    handles.mainGUIhandles.params.gamma = str2num(get(handles.fixedGammaEditText, 'String'));
else
    handles.mainGUIhandles.params.gamma = 1;
end

% Interruptable
if get(handles.interruptableToggle, 'Value') == 1
    handles.mainGUIhandles.params.gui_interrupt = true;
else
    handles.mainGUIhandles.params.gui_interrupt = false;
end

% Verbosity
handles.mainGUIhandles.params.verbose = get(handles.verbosityPopup, 'Value') - 1;

% Update handles
guidata(handles.mainGUI, handles.mainGUIhandles);


% --- Function for making a list of all (compatible) alpha/beta/gamma functions and updating the
%     popup menu
%     func_choise should be 'alpha', 'beta' or 'gamma'
function readFunctions(handles, func_choise)
global DSS_DIRECTORY;
global ALPHA_FUNCTIONS_SORTED;
global BETA_FUNCTIONS_SORTED;
global GAMMA_FUNCTIONS_SORTED;
global CUSTOM_ALPHA_FUNCTIONS;
global CUSTOM_BETA_FUNCTIONS;
global CUSTOM_GAMMA_FUNCTIONS;

switch func_choise
    % Clean popUp menu
    case 'alpha'
        ALPHA_FUNCTIONS_SORTED = [];
        cleanPopup(handles.selectAlphafPopup);
    case 'beta'
        BETA_FUNCTIONS_SORTED = [];
        cleanPopup(handles.selectBetafPopup);
    case 'gamma'
        GAMMA_FUNCTIONS_SORTED = [];
        cleanPopup(handles.selectGammafPopup);
end

dir_struct = dir(DSS_DIRECTORY);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');

k = 1;
found_count = 0;
temp_functions = strcat(func_choise, '_functions');
% Find all the files with '(func_choise)_'-prefix
for i = 1:length(sorted_names)
    if strncmp(sorted_names(i), strcat(func_choise, '_'), length(func_choise)+1)
        handles.(temp_functions){k} = sorted_names(i); 
        k = k + 1;
        found_count = found_count + 1;
    end
end

if found_count == 0
    handles.(temp_functions) = [];
end

if strcmp(func_choise, 'beta')
    supported_betas = betaSupportList(handles);
else
    supported_betas = {};
end

k = 1;
% Make a list of all alpha/beta/gamma functions
for i = 1:length(handles.(temp_functions))
    func_name = char(handles.(temp_functions){i});
    % Remove file extension (.m) if it exists (it should)
    end_char = length(func_name);
    if strcmp(func_name(end_char), '~')
        func_name = '';
    elseif strcmp(func_name(end_char), 'm') & strcmp(func_name(end_char-1), '.')
        func_name = func_name(1:end_char-2);
    end
    if ~isempty(func_name)
        func = str2func(func_name);
    else
        func = [];
    end
    
    beta_support = false;
    % Check denoising support for beta functions
    if ~isempty(supported_betas)
        for j = 1:length(supported_betas)
            if strcmp(func_name, supported_betas{j})
                beta_support = true;
                break;
            end
        end
        if ~beta_support
            func_name = '';
        end
    end
    
    % Try to query function for its parameters
    param_struct = [];
    answered = false;
    try
        param_struct = func([]);
        if ~isempty(func_name)
            answered = true;
        end
    catch
        if ~isempty(func_name)
            % No response...so make it available (responsibility --> user)
            k = addItem(func_choise, func_name, k);
        end
    end
    % If function answered, check its parameters
    if answered
        if ~isempty(param_struct)
            if isfield(param_struct, 'approach')
                numberOfApproaches = length(param_struct.approach);
                for i = 1:numberOfApproaches
                    if strcmp(char(param_struct.approach{i}), handles.mainGUIhandles.params.algorithm)
                        % Match! Make function available..
                        if isfield(param_struct, 'name') 
                            if length(param_struct.name) > 0
                                switch func_choise
                                    case 'alpha'
                                        itemToPopup(handles.selectAlphafPopup, param_struct.name);
                                        ALPHA_FUNCTIONS_SORTED{k} = func_name;
                                        k = k + 1;
                                    case 'beta'
                                        itemToPopup(handles.selectBetafPopup, param_struct.name);
                                        BETA_FUNCTIONS_SORTED{k} = func_name;
                                        k = k + 1;
                                    case 'gamma'
                                        itemToPopup(handles.selectGammafPopup, param_struct.name);
                                        GAMMA_FUNCTIONS_SORTED{k} = func_name;
                                        k = k + 1;
                                end
                            else
                                k = addItem(func_choise, func_name, k);
                            end
                        else
                            k = addItem(func_choise, func_name, k);
                        end
                    end
                end
            end
        else
            if ~isempty(func_name)
                k = addItem(func_choise, func_name, k);
            end
        end
    end
    
end

% Add possible non-standard custom functions to popUp-list
switch func_choise
    case 'alpha'
        n = length(CUSTOM_ALPHA_FUNCTIONS);
        for i = 1:n
            ALPHA_FUNCTIONS_SORTED{k} = CUSTOM_ALPHA_FUNCTIONS{i};
            k = k + 1;
            itemToPopup(handles.selectAlphafPopup, CUSTOM_ALPHA_FUNCTIONS{i});
        end
    case 'beta'
        n = length(CUSTOM_BETA_FUNCTIONS);
        for i = 1:n
            BETA_FUNCTIONS_SORTED{k} = CUSTOM_BETA_FUNCTIONS{i};
            k = k + 1;
            itemToPopup(handles.selectBetafPopup, CUSTOM_BETA_FUNCTIONS{i});
        end
    case 'gamma'
        n = length(CUSTOM_GAMMA_FUNCTIONS);
        for i = 1:n
            GAMMA_FUNCTIONS_SORTED{k} = CUSTOM_GAMMA_FUNCTIONS{i};
            k = k + 1;
            itemToPopup(handles.selectGammafPopup, CUSTOM_GAMMA_FUNCTIONS{i});
        end
end


% --- Support function for readFunctions
%     Adds item to popUp and updates the list of functions
function new_index = addItem(func_choise, func_name, index)

global ALPHA_FUNCTIONS_SORTED;
global BETA_FUNCTIONS_SORTED;
global GAMMA_FUNCTIONS_SORTED;

switch func_choise
    case 'alpha'
        itemToPopup(handles.selectAlphafPopup, func_name);
        ALPHA_FUNCTIONS_SORTED{index} = func_name;
    case 'beta'
        itemToPopup(handles.selectBetafPopup, func_name);
        BETA_FUNCTIONS_SORTED{index} = func_name;
    case 'gamma'
        itemToPopup(handles.selectGammafPopup, func_name);
        GAMMA_FUNCTIONS_SORTED{index} = func_name;
end
new_index = index + 1;


% --- Function for emptying a popUp menu
function cleanPopup(hObject)
% hObject   popUp menu to be emptied
% handles   handles to mainGUI objects
strings = {};
strings{1, 1} = 'Fixed value';
strings{2, 1} = 'Custom';
set(hObject, 'String', strings, 'Value', 1);


% --- Function for adding items to a popUp menu
function itemToPopup(hObject, item)
% hObject   popUp menu to receive item
% item      string to include in popUp
if ~isempty(item)
    strings = get(hObject, 'String');
    numberOfStrings = length(strings);

    % Overwrite 'Custom' with item and add new 'Custom' as the last item
    strings{numberOfStrings, 1} = item;
    strings{numberOfStrings+1, 1} = 'Custom';
    set(hObject, 'String', strings, 'Value', numberOfStrings);
end


% Function for asking what beta functions the current denoising function
% supports
function supported_functions = betaSupportList(handles)

supported_functions = {};
func = handles.mainGUIhandles.params.denf.h;
param_struct = struct();
answered = false;
try
    param_struct = func([], [1 2 3]);
catch
    % No specified support, show all
end

if isfield(param_struct, 'beta')
    if ~isempty(param_struct.beta)
        supported_functions = param_struct.beta;
    end
end