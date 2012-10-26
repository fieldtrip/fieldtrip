function varargout = dss_gui(varargin)
% DSS_GUI M-file for dss_gui.fig
%      DSS_GUI, by itself, creates a new DSS_GUI or raises the existing
%      singleton*.
%
%      H = DSS_GUI returns the handle to a new DSS_GUI or the handle to
%      the existing singleton*.
%
%      DSS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSS_GUI.M with the given input arguments.
%
%      DSS_GUI('Property','Value',...) creates a new DSS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dss_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dss_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      This is the main file for the DSS GUI. Includes functionality of the 
%      main GUI window and calls to other GUI files. Type 'dss_gui' to start.
%      Other DSS GUI files should not be called outside this file.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% Edit the above text to modify the response to help dss_gui

% Last Modified by GUIDE v2.5 28-Feb-2005 10:32:40

% $Id: dss_gui.m,v 1.17 2005/12/07 11:58:28 jaakkos Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dss_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @dss_gui_OutputFcn, ...
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


% --- Executes just before dss_gui is made visible.
function dss_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dss_gui (see VARARGIN)

% Choose default command line output for dss_gui
handles.output = hObject;

% Global variables
global INPUTSIGNAL;
global INPUTSIGNAL_STRING;
INPUTSIGNAL = [];
INPUTSIGNAL_STRING = '';
global DENOISING_FUNCTIONS_SORTED;
DENOISING_FUNCTIONS_SORTED = [];
global PREPRO_FUNCTIONS_SORTED;
PREPRO_FUNCTIONS_SORTED = [];
global ORTHO_FUNCTIONS_SORTED;
ORTHO_FUNCTIONS_SORTED = [];
global LAST_DENOISING_FUNC;
LAST_DENOISING_FUNC = 1;
global LAST_WHITENING_FUNC;
LAST_WHITENING_FUNC = 1;
global LAST_ORTHO_FUNC;
LAST_ORTHO_FUNC = 1;
global DSS_DIRECTORY;
global ALPHA_FUNCTIONS_SORTED;
global BETA_FUNCTIONS_SORTED;
global GAMMA_FUNCTIONS_SORTED;
ALPHA_FUNCTIONS_SORTED = [];
BETA_FUNCTIONS_SORTED = [];
GAMMA_FUNCTIONS_SORTED = [];
global CUSTOM_ALPHA_FUNCTIONS;
global CUSTOM_BETA_FUNCTIONS;
global CUSTOM_GAMMA_FUNCTIONS;
CUSTOM_ALPHA_FUNCTIONS = [];
CUSTOM_BETA_FUNCTIONS = [];
CUSTOM_GAMMA_FUNCTIONS = [];

% -- default functions
def_denoise=struct('h', @denoise_tanh, 'params', '');
def_stop   =struct('h', @default_stop, ...
    'params', struct('epsilon', 1e-4, 'maxiters', 1000));
def_whiten =struct('h', @pre_sphere);
def_ortho  =struct('h', @ortho_default,  'params', '');

% -- parameter & state variables
% Strcture definition: name, description, [default value], [allowed values]
% Not actually used for anything
handles.params_def = {
  {'verbose',      'Output verbosity', 1, {0,1,2,3}},
  {'algorithm',    'DSS algorithm type', 'defl', {'defl', 'symm', 'pca'}},
  {'preprocf',     'Sphering function', def_whiten},
  {'orthof',       'Orthogonalization function', def_ortho},
  {'denf',         'Denoising function', def_denoise},
  {'stopf',        'Stopping criterion function', def_stop},
  {'reportf',      'Vector of iteration reporting functions'},
  {'sdim',         'Projection dimension'},
  {'wdim',         'Dimensions after reduction'},
  {'report',       'Reporting data'},
  {'gui_interrupt','Enable interruption GUI', true},

  {'alphaf',      'Denoise normalization function'},
  {'betaf',       'Denoise spectral shift function'},
  {'gammaf',      'Adaptive gain function'},
  {'adapt_alpha', 'Adaptive (per iteration) alpha', false, {true, false}},
  {'adapt_beta',  'Adaptive (per iteration) beta', false, {true, false}},
  {'adapt_gamma', 'Adaptive (per iteration) gamma', false, {true, false}},
  {'alpha',       'Fixed value for alpha if alphaf is not defined,',1},
  {'beta',        'Fixed value for beta if betaf is not defined,',0},
  {'gamma',       'Fixed value for gamma if gammaf is not defined,',1},
};

% Assigning some default values to params struct
handles.params = struct(...
    'alpha', 1, ...
    'beta', 0, ...
    'gamma', 1, ...
    'adapt_alpha', false, ...
    'adapt_beta', false, ...
    'adapt_gamma', false, ...
    'gui_interrupt', true, ...
    'verbose', 0, ...
    'stopf', def_stop, ...
    'algorithm', 'defl', ...
    'denf', def_denoise, ...
    'orthof', def_ortho, ...
    'preprocf', def_whiten ...
);

% Flag for preprocessing panel default values (true if defaults are selected)
handles.preproDefaultsFlag = true;
% Flag for parameters panel default values (true if defaults are selected)
handles.paramDefaultsFlag = true;

% Name window
set(hObject, 'Name', 'Denoising Source Separation');

% Add some directories to path
dss_dir = which('denss');
% Remove 'denss.m' from the end of the path.
if ~isempty(dss_dir)
    dss_dir = dss_dir(1:(length(dss_dir)-7));
else
    dss_dir = input('Enter full path of the location of the DSS package: ');
end
if ispc
%    misc_dir = strcat(dss_dir, '\misc');
    src_dir = strcat(dss_dir, '\src');
elseif isunix
%    misc_dir = strcat(dss_dir, '/misc');
    src_dir = strcat(dss_dir, '/src');
end
addpath(dss_dir);
% addpath(misc_dir);
addpath(src_dir);
addpath(pwd);
DSS_DIRECTORY = dss_dir;

% Populate denoising popUp with compatible functions
% Compatible here means compatible with Deflation, which is default
readDenoiseFunctions(handles, 'defl');

% Populate whitening popUp with functions
readWhiteningFunctions(handles);

% Populate orthogonalization popUp with functions
readOrthoFunctions(handles);

% Arguments
handles.state_given = false;
if ~isempty(varargin)
    if length(varargin) > 2
        errordlg('Too many input arguments', 'Input Error', 'modal');
        close(hObject);
    end
    for i = 1:length(varargin)
        % Parse arguments
        if isstruct(varargin{i})
            if ~handles.state_given & isempty(INPUTSIGNAL)
                handles.state = varargin{1};
                handles.state_given = true;
                % Enable some elements
                set(handles.plotButton, 'Enable', 'on');
                set(handles.transposeButton, 'Enable', 'on');
                set(handles.startDssButton, 'Enable', 'on');
                set(handles.plotWhitenedButton, 'Enable', 'on');
                guidata(hObject, handles);
                % Update GUI to represent the given state
                INPUTSIGNAL_STRING = 'argument';
                handles = readState(handles, true, true);
            else
                errordlg('Multiple state/input signal arguments. Using first one.', 'Input Warning', 'modal');
            end
        end
        if isnumeric(varargin{i})
            if ~handles.state_given & isempty(INPUTSIGNAL)
                INPUTSIGNAL = varargin{i};
                INPUTSIGNAL_STRING = 'argument';
                % Enable some elements
                set(handles.plotButton, 'Enable', 'on');
                set(handles.transposeButton, 'Enable', 'on');
                set(handles.startDssButton, 'Enable', 'on');
                set(handles.plotWhitenedButton, 'Enable', 'on');
                % Update GUI to represent the given signal
                refreshCounters(handles);
            else
                errordlg('Multiple state/input signal arguments. Using first one.', 'Input Warning', 'modal');
            end
        end
        if ischar(varargin{i})
            % Pro mode, panels 2 and 3 are visible at the start
            if strcmp(varargin{i}, 'pro')
                makeVisible(handles, true, true);
            end
        end
    end
end

% Set default denoise selection in GUI 
handles = setDefaultDenoising(handles);

% Set default preprocessing selection in GUI 
setDefaultPreprocessing(handles);

% Set default orthogonalization selection in GUI
setDefaultOrtho(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dss_gui wait for user response (see UIRESUME)
% uiwait(handles.mainGUI);


% --- Outputs from this function are returned to the command line.
function varargout = dss_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% -------------------------------------------------------
% ------ Start, Exit, Results, Reset all, Help, About ---
% -------------------------------------------------------



% --- Executes on button press in startDssButton.
function startDssButton_Callback(hObject, eventdata, handles)
% hObject    handle to startDssButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Starts DSS processing
% INPUTSIGNAL holds the data we want to process
% handles.params holds current parameters
% When processing is done, reporting function is called

global INPUTSIGNAL;
global DENOISING_FUNCTIONS_SORTED;
global PREPRO_FUNCTIONS_SORTED;
global ORTHO_FUNCTIONS_SORTED;

if isempty(INPUTSIGNAL) & ~handles.state_given
    % No input signal and no state --> raise error
    errordlg('No input signal loaded.', 'Input error', 'modal');
else
    % There is input signal or state
    
    % Gather all the parameters to params-struct and pass it to dss
        if isempty(INPUTSIGNAL)
            source_sig_dim = size(handles.state.X, 1);
        else
            source_sig_dim = size(INPUTSIGNAL, 1);
        end
        
        % WDIM
        if ~isempty(get(handles.wdimInputBox, 'String'))
            handles.params.wdim = str2num(get(handles.wdimInputBox, 'String'));
        else
            handles.params.wdim = handles.params.sdim;
        end
        % wdim can't be bigger than source_sig_dim
        if handles.params.wdim > source_sig_dim | handles.params.wdim < 1
            errordlg('Reduced dimensions must be in range: {1 : source sig dimensions}.', 'Input Error', 'modal');
            return;
        end
        % making sure wdim is integer (rounding if not)
        handles.params.wdim = round(handles.params.wdim);
        guidata(hObject, handles);
        % update rounded value to GUI
        set(handles.wdimInputBox, 'String', num2str(handles.params.wdim));
        
        % SDIM
        if ~isempty(get(handles.sdimEditText, 'String'))
            handles.params.sdim = str2num(get(handles.sdimEditText, 'String'));
        else
            handles.params.sdim = source_sig_dim;
        end
        % sdim can't be bigger than wdim
        if handles.params.sdim > handles.params.wdim | handles.params.sdim < 1
            errordlg('Number of ICs must be in range: {1 : reduced dimensions}.', 'Input Error', 'modal');
            return;
        end
        % making sure sdim is integer (rounding if not)
        handles.params.sdim = round(handles.params.sdim);
        guidata(hObject, handles);
        % update rounded value to GUI
        set(handles.sdimEditText, 'String', num2str(handles.params.sdim));
        
        % WHITENING      
        whiteningfunctions = get(handles.selectWhiteningPopup, 'String');
        index_selected = get(handles.selectWhiteningPopup, 'Value');
        switch whiteningfunctions{index_selected(1)}
            case 'Default whitening'
                handles.params.preprocf.h = @pre_sphere;
            case 'Custom'
                errordlg('Bad selection in whitening menu.', 'Input Error', 'modal');
                return;
            otherwise
                % make a function pointer from the list of function names              
                func_name = char(PREPRO_FUNCTIONS_SORTED{index_selected(1)});
                func = str2func(func_name);
                handles.params.preprocf.h = func;
        end
        
        % APPROACH
        algorithms = get(handles.selectAlgorithmPopup, 'String');
        index_selected = get(handles.selectAlgorithmPopup, 'Value');
        switch algorithms{index_selected(1)}
            case 'Deflation'
                handles.params.algorithm = 'defl';
            case 'Symmetric'
                handles.params.algorithm = 'symm';
            case 'PCA'
                handles.params.algorithm = 'pca';
        end
        
        % DENOISING
        denoisers = get(handles.selectDenoisingPopup, 'String');
        index_selected = get(handles.selectDenoisingPopup, 'Value');
        switch denoisers{index_selected(1)}
            case 'FastICA tanh nonlinearity'
                handles.params.denf.h = @denoise_fica_tanh;
            case 'FastICA kurtosis nonlinearity'
                handles.params.denf.h = @denoise_fica_kurtosis;
            case 'FastICA gaussian nonlinearity'
                handles.params.denf.h = @denoise_fica_gauss;
            case 'FastICA skew nonlinearity'
                handles.params.denf.h = @denoise_fica_skew;
            case 'Supergaussian'
                handles.params.denf.h = @denoise_tanh;
            case 'Smooth tanh'
                handles.params.denf.h = @denoise_smooth_tanh;
            case 'Quasiperiodic averaging'
                handles.params.denf.h = @denoise_avg;
            case 'Mask denoising'
                handles.params.denf.h = @denoise_mask;
            case 'Energy based denoising'
                handles.params.denf.h = @denoise_energy;
            case 'Kurtosis based denoising'
                handles.params.denf.h = @denoise_pow3;
            case 'DCT filter'
                handles.params.denf.h = @denoise_filter;
            case 'Generic filter'
                handles.params.denf.h = @denoise_dct;
            case 'Custom'
                errordlg('Bad selection in denoising menu.', 'Input Error', 'modal');
                return;
            case 'Nothing available'
                errordlg('Bad selection in denoising menu.', 'Input Error', 'modal');
                return;
            otherwise
                % make a function pointer from the list of function names              
                func_name = char(DENOISING_FUNCTIONS_SORTED{index_selected(1)});
                func = str2func(func_name);
                handles.params.denf.h = func;
        end
        
        % ORTHOGONALIZATION
        orthof = get(handles.selectOrthoPopup, 'String');
        index_selected = get(handles.selectOrthoPopup, 'Value');
        switch orthof{index_selected(1)}
            case 'Default orthogonalization'
                handles.params.orthof.h = @ortho_default;
            case 'Quasi-orthogonalization'
                handles.params.orthof.h = @ortho_quasi;
            case 'Custom'
                errordlg('Bad selection in whitening menu.', 'Input Error', 'modal');
                return;
            otherwise
                % make a function pointer from the list of function names              
                func_name = char(ORTHO_FUNCTIONS_SORTED{index_selected(1)});
                func = str2func(func_name);
                handles.params.orthof.h = func;
        end
        
        % STOPPING
        % maximum iterations
        if ~isempty(get(handles.maxItersEditText, 'String'))
            max_iters = str2num(get(handles.maxItersEditText, 'String'));
            max_iters = round(max_iters);
            if max_iters < 1
                errordlg('Max iterations must be at least 1.', 'Input Error', 'modal');
                return;
            end
            handles.params.stopf.params.maxiters = max_iters;
        else
            set(handles.maxItersEditText, 'String', '1000');
            max_iters = str2num(get(handles.maxItersEditText, 'String'));
            handles.params.stopf.params.maxiters = max_iters;
        end
        % epsilon
        if ~isempty(get(handles.epsilonEditText, 'String'))
            epsilon = str2num(get(handles.epsilonEditText, 'String'));
            handles.params.stopf.params.epsilon = epsilon;
        else
            set(handles.epsilonEditText, 'String', '0.1');
            epsilon = str2num(get(handles.epsilonEditText, 'String'));
            handles.params.stopf.params.epsilon = epsilon;
        end
        
        guidata(hObject, handles);

    % DSS call
    if handles.state_given & isempty(INPUTSIGNAL)
        % There is only state
        try
            [state, W, A] = denss(handles.state, handles.params);
            % Save state to handles
            handles.state = state;
            handles.A = A;
            handles.W = W;
        catch
            errordlg(lasterr, 'DSS Error', 'modal');
            return;
        end
    else
        % There is input signal, so use it
        try
            [state, W, A] = denss(INPUTSIGNAL, handles.params);
            % Save state to handles
            handles.state = state;
            handles.A = A;
            handles.W = W;
        catch
            errordlg(lasterr, 'DSS Error', 'modal');
            return;
        end
    end
    
    % Use the unmixing matrix W to the input signal to get independent
    % components
    try
        if handles.state_given & isempty(INPUTSIGNAL)
            S = W * handles.state.X;
        else
            S = W * INPUTSIGNAL;
        end
        handles.S = S;
        guidata(hObject, handles);
        % Enable results-button and show the time of results
        set(handles.resultsButton, 'Enable', 'on');
        time = fix(clock);
        set(handles.resultText, 'Visible', 'on');
        timeInfoText = ['Current results made at ' num2str(time(4)) ':' num2str(time(5)) ':' num2str(time(6))];
        set(handles.resultText, 'String', timeInfoText);
    catch
        errordlg(lasterr, 'DSS Error in using the unmixing matrix', 'modal');
        return;
    end
    
    % Reporting
    fig = figure; clf;
    set(fig, 'Name', 'DSS Results');
    try
        dss_gui_report_result(S, state);
    catch
        errordlg(lasterr, 'DSS Report Function Error', 'modal');
        return;
    end
end


% --- Executes on button press in exitButton.
function exitButton_Callback(hObject, eventdata, handles)
% hObject    handle to exitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clear global variables
clear INPUTSIGNAL;
clear INPUTSIGNAL_STRING;
clear DENOISING_FUNCTIONS_SORTED;
clear PREPRO_FUNCTIONS_SORTED;
clear ORTHO_FUNCTIONS_SORTED;
clear LAST_DENOISING_FUNC;
clear LAST_WHITENING_FUNC;
clear LAST_ORTHO_FUNC;
clear DSS_DIRECTORY;
clear ALPHA_FUNCTIONS_SORTED;
clear BETA_FUNCTIONS_SORTED;
clear GAMMA_FUNCTIONS_SORTED;
clear CUSTOM_ALPHA_FUNCTIONS;
clear CUSTOM_BETA_FUNCTIONS;
clear CUSTOM_GAMMA_FUNCTIONS;

% Close loadGUI if it is open
h = findobj('Tag', 'loadGUI');
if ~isempty(h)
    close(h);
end

% Close advOptGUI if it is open
h = findobj('Tag', 'advOptGUI');
if ~isempty(h)
    close(h);
end

% Close insertfGUI if it is open
h = findobj('Tag', 'insertfGUI');
if ~isempty(h)
    close(h);
end

% Close saveGUI if it is open
h = findobj('Tag', 'saveGUI');
if ~isempty(h)
    close(h);
end

% Close paramsGUI if it is open
h = findobj('Tag', 'paramsGUI');
if ~isempty(h)
    close(h);
end

% Close aboutGUI if it is open
h = findobj('Tag', 'aboutGUI');
if ~isempty(h)
    close(h);
end

% Close mainGUI
close(handles.mainGUI);


% --- Executes on button press in resultsButton.
function resultsButton_Callback(hObject, eventdata, handles)
% hObject    handle to resultsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dss_gui_save(handles);


% --- Executes on button press in resetAllButton.
function resetAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable self
set(hObject, 'Enable', 'off');
% Reset all (resetting to state defaults if state is loaded)
preproResetButton_Callback(handles.preproResetButton, [], handles);
paramResetButton_Callback(handles.paramResetButton, [], handles);


% --- Executes on button press in systemResetButton.
function systemResetButton_Callback(hObject, eventdata, handles)
% hObject    handle to systemResetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable self
set(hObject, 'Enable', 'off');
% Reset all (resetting to system defaults in any case)
handles = setDefaults(handles);
% If advOptions is open, close and re-open it to keep it up to date
h = findobj('Tag', 'advOptGUI');
if ~isempty(h)
    close(h);
    dss_gui_advOptions(handles.mainGUI, handles);
end
guidata(hObject, handles);


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dss_gui_help();


% --- Executes on button press in aboutButton.
function aboutButton_Callback(hObject, eventdata, handles)
% hObject    handle to aboutButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dss_gui_about();



% --------------------------------------------------
% ------ UICONTROLS IN PANEL 1 (LOAD SIGNALS) ------
% --------------------------------------------------



function dataInputBox_Callback(hObject, eventdata, handles)
% hObject    handle to dataInputBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataInputBox as text
%        str2double(get(hObject,'String')) returns contents of dataInputBox as a double


% --- Executes during object creation, after setting all properties.
function dataInputBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataInputBox (see GCBO)
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

% Checking what is written in dataInputBox
% Then checking if there is a variable by that name in the workspace
global INPUTSIGNAL;
global INPUTSIGNAL_STRING;

inputString = get(handles.dataInputBox, 'String');
if ~isempty(inputString)
    vars = evalin('base', 'who');
    vars_length = length(vars);
    for i = 1:vars_length
        if strcmp(inputString, vars(i))
            temp = evalin('base', inputString); 
            % Check if the variable is numeric data and non-empty
            if ~isnumeric(temp)
                if isstruct(temp)
                    % Load state
                    handles.state_given = true;
                    handles.state = temp;
                    INPUTSIGNAL_STRING = inputString;
                    % Update GUI to represent the given state
                    guidata(hObject, handles);
                    handles = readState(handles, true, true);
                    % Enable some elements
                    set(handles.plotButton, 'Enable', 'on');
                    set(handles.transposeButton, 'Enable', 'on');
                    set(handles.startDssButton, 'Enable', 'on');
                    set(handles.systemResetButton, 'Enable', 'on');
                    guidata(hObject, handles);
                    return;
                else
                    errordlg('Variable must be numeric!', 'Input error', 'modal');
                    return;
                end
            end
            if isempty(temp)
                errordlg('Variable is empty.', 'Input error', 'modal');
                return;
            end
            INPUTSIGNAL = temp;
            INPUTSIGNAL_STRING = inputString;
            % If there is a previously loaded state dismiss it and set
            % default parameters
            if handles.state_given
                handles.state_given = false;
                preproResetButton_Callback(handles.preproResetButton, [], handles);
                paramResetButton_Callback(handles.paramResetButton, [], handles);
            end
            refreshCounters(handles);
            % Enable some elements
            set(handles.plotButton, 'Enable', 'on');
            set(handles.transposeButton, 'Enable', 'on');
            set(handles.startDssButton, 'Enable', 'on');
            return;
        end
    end
    % If we get here, nothing was found
    errordlg('No such workspace variable!', 'Input Error', 'modal');
else
    % Call browse if there is no input
    dss_gui_browse(handles.browseWorkspaceButton, handles, handles.mainGUI, [], @readState, @setDefaults);
end


% --- Executes on button press in browseWorkspaceButton.
function browseWorkspaceButton_Callback(hObject, eventdata, handles)
% hObject    handle to browseWorkspaceButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Calling loadVariables-GUI here
dss_gui_browse(hObject, handles, handles.mainGUI, [], @readState, @setDefaults);


% --- Executes on button press in transposeButton.
function transposeButton_Callback(hObject, eventdata, handles)
% hObject    handle to transposeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INPUTSIGNAL;
if ~isempty(INPUTSIGNAL)
    INPUTSIGNAL = transpose(INPUTSIGNAL);
    refreshCounters(handles);
elseif isfield(handles.state, 'X')
    if handles.state_given
        if ~isempty(handles.state.X)
            handles.state.X = transpose(handles.state.X);
            guidata(hObject, handles);
            refreshCounters(handles);
        end
    end
end


% --- Executes on button press in browseFilesButton.
function browseFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to browseFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dss_gui_browse(hObject, handles, handles.mainGUI, [], @readState, @setDefaults);


% --- Executes on button press in plotButton.
function plotButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INPUTSIGNAL;
% This function plots the input data

% State is given
if handles.state_given
    if isfield(handles.state, 'X')
        if ~isempty(handles.state.X)
            dim = size(handles.state.X, 1);
            if dim < 21
                fig = figure; clf;
                set(fig, 'Name', 'Input signals');
                plotindex = 1;
                % Plot
                for signal = 1:dim
                    subplot(dim, 1, plotindex);
                    plot(handles.state.X(signal, :));
                    yscale = max(-min(handles.state.X(signal,:)), max(handles.state.X(signal,:)))*1.3;
                    axis([0 size(handles.state.X,2) -yscale yscale]);
                    plotindex = plotindex + 1;
                end
                return;
            else
                response = dss_gui_modalDlg('Title', ['Confirm plotting of ' num2str(dim) ' signals'], ...
                    'String', 'Plotting that many signals may produce poor results. Are you sure?');
                switch response
                    case 'Yes'
                        fig = figure; clf;
                        set(fig, 'Name', 'Input signals');
                        plotindex = 1;
                        % Plot
                        for signal = 1:dim
                            subplot(dim, 1, plotindex);
                            plot(handles.state.X(signal, :));
                            yscale = max(-min(handles.state.X(signal,:)), max(handles.state.X(signal,:)))*1.3;
                            axis([0 size(handles.state.X,2) -yscale yscale]);
                            plotindex = plotindex + 1;
                        end
                    case 'No'
                        return;
                end
            end
        else
            errordlg('Input signal X in given state is empty.', 'Input Error', 'modal');
            return;
        end
    else
        errordlg('Given state does not have input signal X.', 'Input Error', 'modal');
        return;
    end
    return;
end

% No state, plotting input signal
if ~isempty(INPUTSIGNAL)
    dim = size(INPUTSIGNAL, 1);
    % Check if dim has a good value
    if dim < 21
        fig = figure; clf;
        set(fig, 'Name', 'Input signals');
        plotindex = 1;
        % Plot
        for signal = 1:dim
            subplot(dim, 1, plotindex);
            plot(INPUTSIGNAL(signal, :));
            yscale = max(-min(INPUTSIGNAL(signal,:)), max(INPUTSIGNAL(signal,:)))*1.3;
            axis([0 size(INPUTSIGNAL,2) -yscale yscale]);
            plotindex = plotindex + 1;
        end
    else
        response = dss_gui_modalDlg('Title', ['Confirm plotting of ' num2str(dim) ' signals'], ...
            'String', 'Plotting that many signals may produce poor results. Are you sure?');
        switch response
            case 'Yes'
                fig = figure; clf;
                set(fig, 'Name', 'Input signals');
                plotindex = 1;
                % Plot
                for signal = 1:dim
                    subplot(dim, 1, plotindex);
                    plot(INPUTSIGNAL(signal, :));
                    yscale = max(-min(INPUTSIGNAL(signal,:)), max(INPUTSIGNAL(signal,:)))*1.3;
                    axis([0 size(INPUTSIGNAL,2) -yscale yscale]);
                    plotindex = plotindex + 1;
                end
            case 'No'
                return;
        end
    end
else
    errordlg('No data loaded.', 'Input Error', 'modal');
end


% --- Executes on button press in browseDemosButton.
function browseDemosButton_Callback(hObject, eventdata, handles)
% hObject    handle to browseDemosButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dss_gui_browse(hObject, handles, handles.mainGUI, [], @readState, @setDefaults);



% --------------------------------------------------
% ------ UICONTROLS IN PANEL 2 (PREPROCESSING) -----
% --------------------------------------------------



% --- Executes on button press in dimOptionsButton.
function dimOptionsButton_Callback(hObject, eventdata, handles)
% hObject    handle to dimOptionsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.dimOptionsButton, 'String'), 'Show options')
    % Make elements in preprocessing visible and enable them
    makeVisible(handles, true, false);
else
    % Make elements in preprocessing invisible and disable them
    set(handles.text14, 'Visible', 'off');
    set(handles.text8, 'Visible', 'off');
    set(handles.wdimInputBox, 'Visible', 'off', 'Enable', 'off');
    set(handles.selectWhiteningPopup, 'Visible', 'off', 'Enable', 'off');
    set(handles.plotWhitenedButton, 'Visible', 'off', 'Enable', 'off');
    set(handles.dimOptionsButton, 'String', 'Show options');
    set(handles.preproInfoText, 'Visible', 'on');
    if handles.preproDefaultsFlag
        set(handles.preproInfoText, 'String', 'Default');
    else
        set(handles.preproInfoText, 'String', 'Custom');
    end
end


% --- Executes on button press in preproResetButton.
function preproResetButton_Callback(hObject, eventdata, handles)
% hObject    handle to preproResetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INPUTSIGNAL;
handles.preproDefaultsFlag = true;
set(handles.preproInfoText, 'String', 'Default');
% Disable self
set(hObject, 'Enable', 'off');
% Disable resetAllButton if paramResetButton is disabled
if strcmp(get(handles.paramResetButton, 'Enable'), 'off')
    set(handles.resetAllButton, 'Enable', 'off');
end

if ~handles.state_given
    % There is no state loaded
    setDefaultPreprocessing(handles);
    if ~isempty(INPUTSIGNAL)
        set(handles.wdimInputBox, 'String', size(INPUTSIGNAL, 1));
    else
        set(handles.wdimInputBox, 'String', '');
    end
else
    % There is a loaded state, revert to original state's values
    handles = readState(handles, true, false);
end
    
guidata(hObject, handles);


% --- Executes on selection change in selectWhiteningPopup.
function selectWhiteningPopup_Callback(hObject, eventdata, handles)
% hObject    handle to selectWhiteningPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectWhiteningPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectWhiteningPopup
handles.preproDefaultsFlag = false;
set(handles.preproResetButton, 'Enable', 'on');
set(handles.resetAllButton, 'Enable', 'on');
set(handles.systemResetButton, 'Enable', 'on');
guidata(hObject, handles);
whiteningfunctions = get(hObject, 'String');
index_selected = get(hObject, 'Value');
if strcmp(whiteningfunctions{index_selected(1)}, 'Custom')
    dss_gui_insertFunc(hObject, handles);
else
    dss_gui_funcParams(hObject, handles.mainGUI, handles);
end


% --- Executes during object creation, after setting all properties.
function selectWhiteningPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectWhiteningPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in plotWhitenedButton.
function plotWhitenedButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotWhitenedButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INPUTSIGNAL;
global PREPRO_FUNCTIONS_SORTED;
% This function makes a whitened instance of the input signal
% for itself and plots it

wdim = str2num(get(handles.wdimInputBox, 'String'));
if handles.state_given
    try
        source_sig_dim = size(handles.state.X, 1);
    catch
        source_sig_dim = 0;
        handles.state.X = [];
    end
else
    source_sig_dim = size(INPUTSIGNAL, 1);
end
if isempty(wdim)
    errordlg('Empty input.', 'Input Error', 'modal');
    return;
end
% Check if there is a reasonable amount of signals
good_plot = true;
if wdim > 20
    response = dss_gui_modalDlg('Title', ['Confirm plotting of ' num2str(wdim) ' signals'], ...
        'String', 'Plotting that many signals may produce poor results. Are you sure?');
    switch response
        case 'Yes'
            % Do nothing
        case 'No'
            good_plot = false;
    end
end
if good_plot
    % Check if sdim has a good value
    if wdim <= source_sig_dim & wdim >= 1
        % Making sure wdim is integer (rounding if not)
        wdim = round(wdim);
        % Update corrected value in GUI
        set(handles.wdimInputBox, 'String', num2str(wdim));
        % Get selected preprocessing function
        index_selected = get(handles.selectWhiteningPopup, 'Value');
        prepro_func = str2func(PREPRO_FUNCTIONS_SORTED{index_selected});
        if ~isfield(handles.params.preprocf, 'params')
            handles.params.preprocf.params = [];
        end
        % Call preprocessing function
        if handles.state_given
            if isfield(handles.state, 'X')
                if ~isempty(handles.state.X)
                    [params,wX,wM,dwM] = feval(prepro_func, handles.params.preprocf.params, handles.state.X, wdim);
                else
                    errordlg('Empty data in state.', 'Input Error', 'modal');
                end
            else
                errordlg('No input data in state.', 'Input Error', 'modal');
            end
        else
            if ~isempty(INPUTSIGNAL)
                [params,wX,wM,dwM] = feval(prepro_func, handles.params.preprocf.params, INPUTSIGNAL, wdim);
            else
                errordlg('No input signal loaded.', 'Input Error', 'modal');
            end
        end
        % -------------------
        % Save the whitened data to state struct and activate it
        handles.state.Y = wX;
        handles.state.B = wM;
        handles.state.dB = dwM;
        if ~handles.state_given
            handles.state.X = INPUTSIGNAL;
            INPUTSIGNAL = [];
        end
        handles.state_given = true;
        guidata(hObject, handles);
        % -------------------
        wXdim = size(wX, 1);
        fig = figure; clf;
        set(fig, 'Name', 'Preprocessed data');
        plotindex = 1;
        % Plot
        for signal = 1:wXdim
            subplot(wXdim, 1, plotindex);
            plot(wX(signal, :));
            yscale = max(-min(wX(signal,:)), max(wX(signal,:)))*1.3;
            axis([0 size(wX,2) -yscale yscale]);
            plotindex = plotindex + 1;
        end
    else
        errordlg('Reduced dimensions must be in range: {1 : source signal dimensions}.', 'Input Error', 'modal');
    end
end


function wdimInputBox_Callback(hObject, eventdata, handles)
% hObject    handle to wdimInputBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wdimInputBox as text
%        str2double(get(hObject,'String')) returns contents of wdimInputBox as a double
handles.preproDefaultsFlag = false;
set(handles.preproResetButton, 'Enable', 'on');
set(handles.systemResetButton, 'Enable', 'on');
set(handles.resetAllButton, 'Enable', 'on');
if str2double(get(handles.sdimEditText, 'String')) >= str2double(get(hObject, 'String'))
    set(handles.sdimEditText, 'String', num2str(get(hObject, 'String')));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function wdimInputBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wdimInputBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --------------------------------------------------
% ------ UICONTROLS IN PANEL 3 (CHOOSE PARAMETERS) -
% --------------------------------------------------



% --- Executes on selection change in selectDenoisingPopup.
function selectDenoisingPopup_Callback(hObject, eventdata, handles)
% hObject    handle to selectDenoisingPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectDenoisingPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectDenoisingPopup
handles.paramDefaultsFlag = false;
set(handles.paramResetButton, 'Enable', 'on');
set(handles.resetAllButton, 'Enable', 'on');
set(handles.systemResetButton, 'Enable', 'on');
guidata(hObject, handles);
denoisingfunctions = get(hObject, 'String');
index_selected = get(hObject, 'Value');
if strcmp(denoisingfunctions{index_selected(1)}, 'Custom')
    dss_gui_insertFunc(hObject, handles);
else
    dss_gui_funcParams(hObject, handles.mainGUI, handles);
end


% --- Executes during object creation, after setting all properties.
function selectDenoisingPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectDenoisingPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in selectAlgorithmPopup.
function selectAlgorithmPopup_Callback(hObject, eventdata, handles)
% hObject    handle to selectAlgorithmPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectAlgorithmPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectAlgorithmPopup
handles.paramDefaultsFlag = false;
set(handles.paramResetButton, 'Enable', 'on');
set(handles.resetAllButton, 'Enable', 'on');
set(handles.systemResetButton, 'Enable', 'on');

approaches = get(hObject, 'String');
index_selected = get(hObject, 'Value');
switch approaches{index_selected(1)}
    case 'Deflation'
        handles.params.algorithm = 'defl';
        readDenoiseFunctions(handles, 'defl');
    case 'Symmetric'
        handles.params.algorithm = 'symm';
        readDenoiseFunctions(handles, 'symm');
    case 'PCA'
        handles.params.algorithm = 'pca';
        readDenoiseFunctions(handles, 'pca');
end

% Set default denoise selection in GUI 
handles = setDefaultDenoising(handles);
% Disable alpha, beta and gamma functions if previously one was selected
if isfield(handles.params, 'alphaf')
    handles.params = rmfield(handles.params, 'alphaf');
end
if isfield(handles.params, 'betaf')
    handles.params = rmfield(handles.params, 'betaf');
end
if isfield(handles.params, 'gammaf')
    handles.params = rmfield(handles.params, 'gammaf');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function selectAlgorithmPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectAlgorithmPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in advOptionsButton.
function advOptionsButton_Callback(hObject, eventdata, handles)
% hObject    handle to advOptionsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dss_gui_advOptions(handles.mainGUI, handles);


% --- Executes on selection change in selectOrthoPopup.
function selectOrthoPopup_Callback(hObject, eventdata, handles)
% hObject    handle to selectOrthoPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectOrthoPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectOrthoPopup
handles.paramDefaultsFlag = false;
set(handles.paramResetButton, 'Enable', 'on');
set(handles.resetAllButton, 'Enable', 'on');
set(handles.systemResetButton, 'Enable', 'on');
guidata(hObject, handles);
orthofunctions = get(hObject, 'String');
index_selected = get(hObject, 'Value');
if strcmp(orthofunctions{index_selected(1)}, 'Custom')
    dss_gui_insertFunc(hObject, handles);
else
    dss_gui_funcParams(hObject, handles.mainGUI, handles);
end


% --- Executes during object creation, after setting all properties.
function selectOrthoPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectOrthoPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function maxItersEditText_Callback(hObject, eventdata, handles)
% hObject    handle to maxItersEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxItersEditText as text
%        str2double(get(hObject,'String')) returns contents of maxItersEditText as a double
handles.paramDefaultsFlag = false;
set(handles.paramResetButton, 'Enable', 'on');
set(handles.resetAllButton, 'Enable', 'on');
set(handles.systemResetButton, 'Enable', 'on');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function maxItersEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxItersEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function epsilonEditText_Callback(hObject, eventdata, handles)
% hObject    handle to epsilonEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilonEditText as text
%        str2double(get(hObject,'String')) returns contents of epsilonEditText as a double
handles.paramDefaultsFlag = false;
set(handles.paramResetButton, 'Enable', 'on');
set(handles.resetAllButton, 'Enable', 'on');
set(handles.systemResetButton, 'Enable', 'on');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function epsilonEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilonEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in paramOptionsButton.
function paramOptionsButton_Callback(hObject, eventdata, handles)
% hObject    handle to paramOptionsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.paramOptionsButton, 'String'), 'Show options')
    % Make elements in choose parameters visible and enable them
    makeVisible(handles, false, true);
else
    % Make elements in choose parameters invisible and disable them
    set(handles.uipanel6, 'Visible', 'off');
    set(handles.text6, 'Visible', 'off');
    set(handles.text7, 'Visible', 'off');
    set(handles.sdimEditText, 'Visible', 'off', 'Enable', 'off');
    set(handles.text9, 'Visible', 'off');
    set(handles.text11, 'Visible', 'off');
    set(handles.text12, 'Visible', 'off');
    set(handles.text13, 'Visible', 'off');
    set(handles.selectOrthoPopup, 'Visible', 'off', 'Enable', 'off');
    set(handles.selectAlgorithmPopup, 'Visible', 'off', 'Enable', 'off');
    set(handles.selectDenoisingPopup, 'Visible', 'off', 'Enable', 'off');
    set(handles.maxItersEditText, 'Visible', 'off', 'Enable', 'off');
    set(handles.epsilonEditText, 'Visible', 'off', 'Enable', 'off');
    set(handles.advOptionsButton, 'Visible', 'off', 'Enable', 'off');
    set(handles.paramOptionsButton, 'String', 'Show options');
    set(handles.paramInfoText, 'Visible', 'on');
    if handles.paramDefaultsFlag
        set(handles.paramInfoText, 'String', 'Default');
    else
        set(handles.paramInfoText, 'String', 'Custom');
    end
end


% --- Executes on button press in paramResetButton.
function paramResetButton_Callback(hObject, eventdata, handles)
% hObject    handle to paramResetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INPUTSIGNAL;
handles.paramDefaultsFlag = true;
set(handles.paramInfoText, 'String', 'Default');
% Disable self
set(hObject, 'Enable', 'off');
% Disable resetAllButton if preproResetButton is disabled
if strcmp(get(handles.preproResetButton, 'Enable'), 'off')
    set(handles.resetAllButton, 'Enable', 'off');
end

if ~handles.state_given
    % No state loaded, load default parameters
    
    % Set default denoise selection in GUI 
    handles = setDefaultDenoising(handles);
    
    % Set default orthogonalization selection in GUI 
    setDefaultOrtho(handles);
    
    set(handles.selectAlgorithmPopup, 'Value', 1);
    set(handles.maxItersEditText, 'String', '1000');
    set(handles.epsilonEditText, 'String', '0.1');
    if ~isempty(INPUTSIGNAL)
        set(handles.sdimEditText, 'String', size(INPUTSIGNAL, 1));
    else
        set(handles.sdimEditText, 'String', '');
    end
    % Defaults for adv. options
    handles.params.gui_interrupt = true;
    handles.params.verbose = 0;
    handles.params.alpha = 1;
    handles.params.adapt_alpha = false;
    handles.params.beta = 0;
    handles.params.adapt_beta = false;
    handles.params.gamma = 1;
    handles.params.adapt_gamma = false;
    if isfield(handles.params, 'alphaf')
        handles.params = rmfield(handles.params, 'alphaf');
    end
    if isfield(handles.params, 'betaf')
        handles.params = rmfield(handles.params, 'betaf');
    end
    if isfield(handles.params, 'gammaf')
        handles.params = rmfield(handles.params, 'gammaf');
    end
else
    % State is loaded, load original state parameters
    handles = readState(handles, false, true);
end

% If advOptions is open, close and re-open it to keep it up to date
h = findobj('Tag', 'advOptGUI');
if ~isempty(h)
    close(h);
    dss_gui_advOptions(handles.mainGUI, handles);
end

guidata(hObject, handles);


function sdimEditText_Callback(hObject, eventdata, handles)
% hObject    handle to sdimEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sdimEditText as text
%        str2double(get(hObject,'String')) returns contents of sdimEditText as a double
handles.paramDefaultsFlag = false;
set(handles.paramResetButton, 'Enable', 'on');
set(handles.resetAllButton, 'Enable', 'on');
set(handles.systemResetButton, 'Enable', 'on');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sdimEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sdimEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --------------------------------------------------
% ------ SUPPORT FUNCTIONS -------------------------
% --------------------------------------------------



% --- Goes through all files with 'denoise'-prefix in the current directory
%     and constructs the list of valid denoising functions for current
%     approach (deflation, symmetric or PCA)
function readDenoiseFunctions(handles, approach)

global DENOISING_FUNCTIONS_SORTED;
global DSS_DIRECTORY;
DENOISING_FUNCTIONS_SORTED = [];
dummy_signal = [1 2 3];
% Clean popUp menu
cleanPopup(handles.selectDenoisingPopup);

dir_struct = dir(DSS_DIRECTORY);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');

k = 1;
found_count = 0;
% Find all the files with 'denoise_'-prefix
for i = 1:length(sorted_names)
    if strncmp(sorted_names(i), 'denoise_', 8)
       handles.denoise_functions{k} = sorted_names(i); 
       k = k + 1;
       found_count = found_count + 1;
    end
end

if found_count == 0
    handles.denoise_functions = [];
end

k = 1;
% Query for supported approaches
for i = 1:length(handles.denoise_functions)
    func_name = char(handles.denoise_functions{i});
    % Remove file extension (.m) if it exists (it should)
    end_char = length(func_name);
    % Skip if the file is a copy (.m~)
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
    % Try to query function for its parameters
    param_struct = [];
    answered = false;
    try
        param_struct = func([], dummy_signal);
        if ~isempty(param_struct)
            answered = true;
        end
    catch
        % No response...so make it available (responsibility --> user)
        if ~isempty(func_name)
            itemToPopup(handles.selectDenoisingPopup, func_name);
            DENOISING_FUNCTIONS_SORTED{k} = func_name;
            k = k + 1;
        end
    end
    % If function answered, check its parameters
    if answered
        if ~isempty(param_struct)
            if isfield(param_struct, 'approach')
                numberOfApproaches = length(param_struct.approach);
                for i = 1:numberOfApproaches
                    if strcmp(char(param_struct.approach{i}), approach)
                        % Match! Make function available..
                        if isfield(param_struct, 'name') 
                            if length(param_struct.name) > 0
                                itemToPopup(handles.selectDenoisingPopup, param_struct.name);
                                DENOISING_FUNCTIONS_SORTED{k} = func_name;
                                k = k + 1;
                            else
                                itemToPopup(handles.selectDenoisingPopup, func_name);
                                DENOISING_FUNCTIONS_SORTED{k} = func_name;
                                k = k + 1;
                            end
                        else
                            itemToPopup(handles.selectDenoisingPopup, func_name);
                            DENOISING_FUNCTIONS_SORTED{k} = func_name;
                            k = k + 1;
                        end
                    end
                end
            end
        else
            if ~isempty(func_name)
                itemToPopup(handles.selectDenoisingPopup, func_name);
                DENOISING_FUNCTIONS_SORTED{k} = func_name;
                k = k + 1;
            end
        end
    end
    
end

% Add possible non-standard custom functions to popUp-list 
if isfield(handles, 'customDenoisingFunctions')
    n = length(handles.customDenoisingFunctions);
    for i = 1:n
        DENOISING_FUNCTIONS_SORTED{k} = handles.customDenoisingFunctions{i};
        k = k + 1;
        itemToPopup(handles.selectDenoisingPopup, handles.customDenoisingFunctions{i});
    end
end

% TODO: Sort pop up.....maybe name the files better?
%sorted_popUp = sortrows(get(mainGUIhandles.selectDenoisingPopup, 'String'));
%set(mainGUIhandles.selectDenoisingPopup, 'String', sorted_popUp);


% --- Function for making panels 2 and 3 visible
%     (preprocessing and choose parameters)
function makeVisible(handles, panel_2, panel_3)
% panel_2   Makes panel_2 visible if true
% panel_3   Makes paneö_3 visible if true
if panel_2
    set(handles.text14, 'Visible', 'on');
    set(handles.text8, 'Visible', 'on');
    set(handles.wdimInputBox, 'Visible', 'on', 'Enable', 'on');
    set(handles.selectWhiteningPopup, 'Visible', 'on', 'Enable', 'on');
    set(handles.plotWhitenedButton, 'Visible', 'on', 'Enable', 'on');
    set(handles.preproInfoText, 'Visible', 'off');
    set(handles.dimOptionsButton, 'String', 'Hide options');
end
if panel_3
    set(handles.sdimEditText, 'Visible', 'on', 'Enable', 'on');
    set(handles.uipanel6, 'Visible', 'on');
    set(handles.text6, 'Visible', 'on');
    set(handles.text7, 'Visible', 'on');
    set(handles.text9, 'Visible', 'on');
    set(handles.text11, 'Visible', 'on');
    set(handles.text12, 'Visible', 'on');
    set(handles.text13, 'Visible', 'on');
    set(handles.selectOrthoPopup, 'Visible', 'on', 'Enable', 'on');
    set(handles.selectAlgorithmPopup, 'Visible', 'on', 'Enable', 'on');
    set(handles.selectDenoisingPopup, 'Visible', 'on', 'Enable', 'on');
    set(handles.maxItersEditText, 'Visible', 'on', 'Enable', 'on');
    set(handles.epsilonEditText, 'Visible', 'on', 'Enable', 'on');
    set(handles.advOptionsButton, 'Visible', 'on', 'Enable', 'on');
    set(handles.paramInfoText, 'Visible', 'off');
    set(handles.paramOptionsButton, 'String', 'Hide options'); 
end


% --- Function for emptying a popUp menu
function cleanPopup(hObject)
% hObject   popUp menu to be emptied
% handles   handles to mainGUI objects
strings = {};
strings{1, 1} = 'Nothing available';
strings{2, 1} = 'Custom';
set(hObject, 'String', strings, 'Value', 1);


% --- Function for adding items to a popUp menu
function itemToPopup(hObject, item)
% hObject   popUp menu to receive item
% item      string to include in popUp
strings = get(hObject, 'String');
numberOfStrings = length(strings);
if strcmp(char(strings{1, 1}), 'Nothing available')
    strings{1, 1} = item;
    set(hObject, 'String', strings, 'Value', 1);
else
    % Overwrite 'Custom' with item and add new 'Custom' as the last item
    strings{numberOfStrings, 1} = item;
    strings{numberOfStrings+1, 1} = 'Custom';
    set(hObject, 'String', strings, 'Value', numberOfStrings);
end


% --- Function for constructing list of whitening functions
function readWhiteningFunctions(handles)
global PREPRO_FUNCTIONS_SORTED;
global DSS_DIRECTORY;
PREPRO_FUNCTIONS_SORTED = [];

% Clean popUp menu
cleanPopup(handles.selectWhiteningPopup);

dir_struct = dir(DSS_DIRECTORY);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');

k = 1;
found_count = 0;
% Find all the files with 'pre_'-prefix
for i = 1:length(sorted_names)
    if strncmp(sorted_names(i), 'pre_', 4)
       handles.whitening_functions{k} = sorted_names(i); 
       k = k + 1;
       found_count = found_count + 1;
    end
end

if found_count == 0
    handles.whitening_functions = [];
end

k = 1;
% Make a list of all preprocessing functions
for i = 1:length(handles.whitening_functions)
    func_name = char(handles.whitening_functions{i});
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
            itemToPopup(handles.selectWhiteningPopup, func_name);
            PREPRO_FUNCTIONS_SORTED{k} = func_name;
            k = k + 1;
        end
    end
    % If function answered, check its parameters
    if answered
        if ~isempty(param_struct)
            if isfield(param_struct, 'name') 
                if length(param_struct.name) > 0
                    itemToPopup(handles.selectWhiteningPopup, param_struct.name);
                    PREPRO_FUNCTIONS_SORTED{k} = func_name;
                    k = k + 1;
                else
                    itemToPopup(handles.selectWhiteningPopup, func_name);
                    PREPRO_FUNCTIONS_SORTED{k} = func_name;
                    k = k + 1;
                end
            else
                itemToPopup(handles.selectWhiteningPopup, func_name);
                PREPRO_FUNCTIONS_SORTED{k} = func_name;
                k = k + 1;
            end
        else
            if ~isempty(func_name)
                itemToPopup(handles.selectDenoisingPopup, func_name);
                PREPRO_FUNCTIONS_SORTED{k} = func_name;
                k = k + 1;
            end
        end
    end
    
end


% --- Function for constructing list of orthogonalization functions
function readOrthoFunctions(handles)
global ORTHO_FUNCTIONS_SORTED;
global DSS_DIRECTORY;
ORTHO_FUNCTIONS_SORTED = [];

% Clean popUp menu
cleanPopup(handles.selectOrthoPopup);

dir_struct = dir(DSS_DIRECTORY);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');

k = 1;
found_count = 0;
% Find all the files with 'ortho_'-prefix
for i = 1:length(sorted_names)
    if strncmp(sorted_names(i), 'ortho_', 6)
       handles.ortho_functions{k} = sorted_names(i); 
       k = k + 1;
       found_count = found_count + 1;
    end
end

if found_count == 0
    handles.ortho_functions = [];
end

k = 1;
% Make a list of all ortho functions
for i = 1:length(handles.ortho_functions)
    func_name = char(handles.ortho_functions{i});
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
            itemToPopup(handles.selectOrthoPopup, func_name);
            ORTHO_FUNCTIONS_SORTED{k} = func_name;
            k = k + 1;
        end
    end
    % If function answered, check its parameters
    if answered
        if ~isempty(param_struct)
            if isfield(param_struct, 'name') 
                if length(param_struct.name) > 0
                    itemToPopup(handles.selectOrthoPopup, param_struct.name);
                    ORTHO_FUNCTIONS_SORTED{k} = func_name;
                    k = k + 1;
                else
                    itemToPopup(handles.selectOrthoPopup, func_name);
                    ORTHO_FUNCTIONS_SORTED{k} = func_name;
                    k = k + 1;
                end
            else
                itemToPopup(handles.selectOrthoPopup, func_name);
                ORTHO_FUNCTIONS_SORTED{k} = func_name;
                k = k + 1;
            end
        else
            if ~isempty(func_name)
                itemToPopup(handles.selectOrthoPopup, func_name);
                ORTHO_FUNCTIONS_SORTED{k} = func_name;
                k = k + 1;
            end
        end
    end
    
end


% --- Refresh data in GUI when needed
function refreshCounters(handles)

global INPUTSIGNAL;
global INPUTSIGNAL_STRING;

% Update signal dimension counters and informText...

if ~isempty(INPUTSIGNAL)
    [rows, columns] = size(INPUTSIGNAL);
    set(handles.numberOfSignalsText, 'String', num2str(rows));
    set(handles.numberOfSamplesText, 'String', num2str(columns));
    set(handles.informText, 'String', ['Signal ' INPUTSIGNAL_STRING ' loaded.']);
    set(handles.sdimEditText, 'String', num2str(rows));
    set(handles.wdimInputBox, 'String', num2str(rows));
else
    if isfield(handles, 'state')
        if isfield(handles.state, 'X')
            if ~isempty(handles.state.X)
                if handles.state_given
                    set(handles.informText, 'String', ['State ' INPUTSIGNAL_STRING ' loaded.']);
                    set(handles.numberOfSignalsText, 'String', num2str(size(handles.state.X, 1)));
                    set(handles.numberOfSamplesText, 'String', num2str(size(handles.state.X, 2)));
                    set(handles.wdimInputBox, 'String', num2str(handles.state.wdim));
                    set(handles.sdimEditText, 'String', num2str(handles.state.sdim));
                else
                    set(handles.numberOfSignalsText, 'String', '0');
                    set(handles.numberOfSamplesText, 'String', '0');
                end
            end
        end
    end
end
if isempty(INPUTSIGNAL) & ~handles.state_given
    set(handles.numberOfSignalsText, 'String', '0');
    set(handles.numberOfSamplesText, 'String', '0');
    set(handles.wdimInputBox, 'String', '');
    set(handles.sdimEditText, 'String', '');
end


% --- Function for reading given state parameter and changing selections in
%     the GUI accordingly
function handles = readState(handles, prepro, params)
% prepro    if true, reads parameters that have elements in preprocessing
%           panel, otherwise not
% params    if true, reads parameters that have elements in choose
%           parameters panel, otherwise not
% Returns new handles struct

global PREPRO_FUNCTIONS_SORTED;
global DENOISING_FUNCTIONS_SORTED;
global ORTHO_FUNCTIONS_SORTED;
global INPUTSIGNAL_STRING;
global CUSTOM_ALPHA_FUNCTIONS;
global CUSTOM_BETA_FUNCTIONS;
global CUSTOM_GAMMA_FUNCTIONS;

% Update informText to show that a state is loaded
set(handles.informText, 'String', ['State ' INPUTSIGNAL_STRING ' loaded.']);

% Update dimensions and samples -display in GUI
if isfield(handles.state, 'X')
    if ~isempty(handles.state.X)
        set(handles.numberOfSignalsText, 'String', num2str(size(handles.state.X, 1)));
        set(handles.numberOfSamplesText, 'String', num2str(size(handles.state.X, 2)));
    end
end

% Read all the fields in the state struct and update GUI selections
% accordingly

if prepro
    % WDIM
    if isfield(handles.state, 'wdim')
        handles.params.wdim = handles.state.wdim;
        set(handles.wdimInputBox, 'String', num2str(handles.state.wdim));
    end
    % WHITENING
    if isfield(handles.state, 'preprocf')
        handles.params.preprocf = handles.state.preprocf;
        found = false;
        numberOfWhiteningFunctions = length(PREPRO_FUNCTIONS_SORTED);
        for i = 1:numberOfWhiteningFunctions
            if strcmp(func2str(handles.state.preprocf.h), PREPRO_FUNCTIONS_SORTED{i})
                set(handles.selectWhiteningPopup, 'Value', i);
                found = true;
                break;
            end
        end
        if ~found
            itemToPopup(handles.selectWhiteningPopup, func2str(handles.state.preprocf.h));
        end
    end
end

if params
    % SDIM
    if isfield(handles.state, 'sdim')
        handles.params.sdim = handles.state.sdim;
        set(handles.sdimEditText, 'String', num2str(handles.state.sdim));
    end
    % VERBOSITY
    if isfield(handles.state, 'verbose')
        handles.params.verbose = handles.state.verbose;
    end
    % GUI INTERRUPT
    if isfield(handles.state, 'gui_interrupt')
        handles.params.gui_interrupt = handles.state.gui_interrupt;
    end
    % APPROACH
    if isfield(handles.state, 'algorithm')
        handles.params.algorithm = handles.state.algorithm;
        switch handles.state.algorithm
            case 'defl'
                set(handles.selectAlgorithmPopup, 'Value', 1);
            case 'symm'
                set(handles.selectAlgorithmPopup, 'Value', 2);
            case 'pca'
                set(handles.selectAlgorithmPopup, 'Value', 3);
        end
    end
    % DENOISING
    if isfield(handles.state, 'denf')
        handles.params.denf = handles.state.denf;
        readDenoiseFunctions(handles, handles.params.algorithm);
        found = false;
        numberOfDenoiseFunctions = length(DENOISING_FUNCTIONS_SORTED);
        for i = 1:numberOfDenoiseFunctions
            if strcmp(func2str(handles.state.denf.h), DENOISING_FUNCTIONS_SORTED{i})
                set(handles.selectDenoisingPopup, 'Value', i);
                found = true;
                break;
            end
        end
        if ~found
            itemToPopup(handles.selectDenoisingPopup, func2str(handles.state.denf.h));
            DENOISING_FUNCTIONS_SORTED{numberOfDenoiseFunctions+1} = func2str(handles.state.denf.h);
            if isfield(handles, 'customDenoisingFunctions')
                n = length(handles.customDenoisingFunctions);
                handles.customDenoisingFunctions{n+1} = func2str(handles.state.denf.h);
            else
                handles.customDenoisingFunctions = {};
                handles.customDenoisingFunctions{1} = func2str(handles.state.denf.h);
            end
        end
    end
    % ORTHOGONALIZATION
    if isfield(handles.state, 'orthof')
        handles.params.orthof = handles.state.orthof;
        found = false;
        numberOfOrthoFunctions = length(ORTHO_FUNCTIONS_SORTED);
        for i = 1:numberOfOrthoFunctions
            if strcmp(func2str(handles.state.orthof.h), ORTHO_FUNCTIONS_SORTED{i})
                set(handles.selectOrthoPopup, 'Value', i);
                found = true;
                break;
            end
        end
        if ~found
            itemToPopup(handles.selectOrthoPopup, func2str(handles.state.orthof.h));
        end
    end
    % STOPPING
    if isfield(handles.state, 'stopf')
        if isfield(handles.state.stopf, 'params')
            if isfield(handles.state.stopf.params, 'epsilon')
                handles.params.stopf.params.epsilon = handles.state.stopf.params.epsilon;
                set(handles.epsilonEditText, 'String', num2str(handles.state.stopf.params.epsilon));
            end
            if isfield(handles.state.stopf.params, 'maxiters')
                handles.params.stopf.params.maxiters = handles.state.stopf.params.maxiters;
                set(handles.maxItersEditText, 'String', num2str(handles.state.stopf.params.maxiters));
            end
        end
    end
    % ALPHA, BETA, GAMMA
    % Functions
    found = false;
    if isfield(handles.state, 'alphaf')
        if ~isempty(handles.state.alphaf)
            handles.params.alphaf = handles.state.alphaf;
            for i = 1:length(CUSTOM_ALPHA_FUNCTIONS)
                if strcmp(CUSTOM_ALPHA_FUNCTIONS{i}, func2str(handles.state.alphaf.h))
                    found = true;
                    break;
                end
            end
            if ~found
                CUSTOM_ALPHA_FUNCTIONS{length(CUSTOM_ALPHA_FUNCTIONS)+1} = func2str(handles.state.alphaf.h);
            end
        end
    end
    found = false;
    if isfield(handles.state, 'betaf')
        if ~isempty(handles.state.betaf)
            handles.params.betaf = handles.state.betaf;
            for i = 1:length(CUSTOM_BETA_FUNCTIONS)
                if strcmp(CUSTOM_BETA_FUNCTIONS{i}, func2str(handles.state.betaf.h))
                    found = true;
                    break;
                end
            end
            if ~found
                CUSTOM_BETA_FUNCTIONS{length(CUSTOM_BETA_FUNCTIONS)+1} = func2str(handles.state.betaf.h);
            end
        end
    end
    found = false;
    if isfield(handles.state, 'gammaf')
        if ~isempty(handles.state.gammaf)
            handles.params.gammaf = handles.state.gammaf;
            for i = 1:length(CUSTOM_GAMMA_FUNCTIONS)
                if strcmp(CUSTOM_GAMMA_FUNCTIONS{i}, func2str(handles.state.gammaf.h))
                    found = true;
                    break;
                end
            end
            if ~found
                CUSTOM_GAMMA_FUNCTIONS{length(CUSTOM_GAMMA_FUNCTIONS)+1} = func2str(handles.state.gammaf.h);
            end
        end
    end
    % Fixed values
    if isfield(handles.state, 'alpha')
        if ~isempty(handles.state.alpha)
            handles.params.alpha = handles.state.alpha;
        end
    end
    if isfield(handles.state, 'beta')
        if ~isempty(handles.state.beta)
            handles.params.beta = handles.state.beta;
        end
    end
    if isfield(handles.state, 'gamma')
        if ~isempty(handles.state.gamma)
            handles.params.gamma = handles.state.gamma;
        end
    end
    % Adaptives
    if isfield(handles.state, 'adapt_alpha')
        if ~isempty(handles.state.adapt_alpha)
            handles.params.adapt_alpha = handles.state.adapt_alpha;
        end
    end
    if isfield(handles.state, 'adapt_beta')
        if ~isempty(handles.state.adapt_beta)
            handles.params.adapt_beta = handles.state.adapt_beta;
        end
    end
    if isfield(handles.state, 'adapt_gamma')
        if ~isempty(handles.state.adapt_gamma)
            handles.params.adapt_gamma = handles.state.adapt_gamma;
        end
    end
end


% --- Function for setting default values for all GUI selections
function handles = setDefaults(handles)
global INPUTSIGNAL;

    setDefaultPreprocessing(handles);
    
    if ~handles.state_given
        if ~isempty(INPUTSIGNAL)
            set(handles.wdimInputBox, 'String', size(INPUTSIGNAL, 1));
        else
            set(handles.wdimInputBox, 'String', '0');
        end
    else
        if isfield(handles.state, 'X')
            if ~isempty(handles.state.X)
                set(handles.wdimInputBox, 'String', size(handles.state.X, 1));
            else
                set(handles.wdimInputBox, 'String', '0');
            end
        else
            set(handles.wdimInputBox, 'String', '0');
        end
    end

    handles = setDefaultDenoising(handles);
    setDefaultOrtho(handles);
    
    set(handles.selectAlgorithmPopup, 'Value', 1);
    set(handles.maxItersEditText, 'String', '1000');
    set(handles.epsilonEditText, 'String', '0.1');
    if ~handles.state_given
        if ~isempty(INPUTSIGNAL)
            set(handles.sdimEditText, 'String', size(INPUTSIGNAL, 1));
        else
            set(handles.sdimEditText, 'String', '0');
        end
    else
        if isfield(handles.state, 'X')
            if ~isempty(handles.state.X)
                set(handles.sdimEditText, 'String', size(handles.state.X, 1));
            else
                set(handles.sdimEditText, 'String', '0');
            end
        else
            set(handles.sdimEditText, 'String', '0');
        end
    end
    % Defaults for adv. options
    handles.params.gui_interrupt = true;
    handles.params.verbose = 0;
    handles.params.alpha = 1;
    handles.params.adapt_alpha = false;
    handles.params.beta = 0;
    handles.params.adapt_beta = false;
    handles.params.gamma = 1;
    handles.params.adapt_gamma = false;
    if isfield(handles.params, 'alphaf')
        handles.params = rmfield(handles.params, 'alphaf');
    end
    if isfield(handles.params, 'betaf')
        handles.params = rmfield(handles.params, 'betaf');
    end
    if isfield(handles.params, 'gammaf')
        handles.params = rmfield(handles.params, 'gammaf');
    end
    

% --- Function for setting default value in preprocessing popUp    
function setDefaultPreprocessing(handles)
global LAST_WHITENING_FUNC;

% Set default preprocessing selection in GUI 

    % Select Default whitening as default
    prepro_functions = get(handles.selectWhiteningPopup, 'String');
    for i = 1:length(prepro_functions)
        if strcmp(prepro_functions{i}, 'Default whitening')
            set(handles.selectWhiteningPopup, 'Value', i);
            LAST_WHITENING_FUNC = i;
            break;
        end
    end


% --- Function for setting default value in denoising popUp    
function handles = setDefaultDenoising(handles)
global LAST_DENOISING_FUNC;

% Set default denoise selection in GUI 

    % Select Supergaussian as default
    denoising_functions = get(handles.selectDenoisingPopup, 'String');
    for i = 1:length(denoising_functions)
        if strcmp(denoising_functions{i}, 'Supergaussian')
            set(handles.selectDenoisingPopup, 'Value', i);
            LAST_DENOISING_FUNC = i;
            handles.params.denf.h = @denoise_tanh;
            handles.params.denf.params = [];
            % Select local beta
            handles.params.betaf.h = @beta_tanh;
            guidata(handles.mainGUI, handles);
            break;
        end
    end


% --- Function for setting default value in orthogonalization popUp    
function setDefaultOrtho(handles)
global LAST_ORTHO_FUNC;

% Set default orthogonalization selection in GUI 

    % Select Default orthogonalization as default
    ortho_functions = get(handles.selectOrthoPopup, 'String');
    for i = 1:length(ortho_functions)
        if strcmp(ortho_functions{i}, 'Default orthogonalization')
            set(handles.selectOrthoPopup, 'Value', i);
            LAST_ORTHO_FUNC = i;
            break;
        end
    end