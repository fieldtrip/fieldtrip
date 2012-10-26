function varargout = dss_gui_funcParams(varargin)
% DSS_GUI_FUNCPARAMS M-file for dss_gui_funcParams.fig
%      DSS_GUI_FUNCPARAMS, by itself, creates a new DSS_GUI_FUNCPARAMS or raises the existing
%      singleton*.
%
%      H = DSS_GUI_FUNCPARAMS returns the handle to a new DSS_GUI_FUNCPARAMS or the handle to
%      the existing singleton*.
%
%      DSS_GUI_FUNCPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSS_GUI_FUNCPARAMS.M with the given input arguments.
%
%      DSS_GUI_FUNCPARAMS('Property','Value',...) creates a new DSS_GUI_FUNCPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dss_gui_funcParams_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dss_gui_funcParams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Includes functionality of the function parameters window.
%      This file is used by other DSS files and should not be used alone.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% Edit the above text to modify the response to help dss_gui_funcParams

% Last Modified by GUIDE v2.5 24-Feb-2005 14:52:43

% $Id: dss_gui_funcParams.m,v 1.12 2005/05/16 12:31:58 kosti Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dss_gui_funcParams_OpeningFcn, ...
                   'gui_OutputFcn',  @dss_gui_funcParams_OutputFcn, ...
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


% --- Executes just before dss_gui_funcParams is made visible.
function dss_gui_funcParams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dss_gui_funcParams (see VARARGIN)

% Choose default command line output for dss_gui_funcParams
handles.output = hObject;

global DENOISING_FUNCTIONS_SORTED;
global PREPRO_FUNCTIONS_SORTED;
global ORTHO_FUNCTIONS_SORTED;

% Name dialog
set(hObject, 'Name', 'Function parameters');

% Make dialog modal
set(hObject, 'WindowStyle','modal');

% Handles from mainGUI
handles.caller = varargin{1};
handles.mainGUI = varargin{2};
handles.mainGUIhandles = varargin{3};

% Find out who called this dialog
guidata(hObject, handles);
handles.callerTag = get(handles.caller, 'Tag');
h = findobj('Tag', handles.callerTag);

% Get the value of mainGUI popUp that called
handles.valueOfPopup = get(h, 'Value');

% Get functions name based on who called
switch handles.callerTag 
    case 'selectDenoisingPopup'
        func_name = char(DENOISING_FUNCTIONS_SORTED{handles.valueOfPopup});
        handles.func_abbr = 'denf';
    case 'selectOrthoPopup'
        func_name = char(ORTHO_FUNCTIONS_SORTED{handles.valueOfPopup});
        handles.func_abbr = 'orthof';
    case 'selectWhiteningPopup'
        func_name = char(PREPRO_FUNCTIONS_SORTED{handles.valueOfPopup});
        handles.func_abbr = 'preprocf';
end

set(handles.fnameText, 'String', func_name);
func = str2func(func_name);
handles.func = func;
handles.answered = false;
handles.asked_params = [];
handles.params = [];
handles.numberOfParams = 0;
handles.numberOfValues = 0;
handles.varnames = {};
% Clean up
cleanPopup(handles.parameterPopup);
% Find out what parameters the function has
try
    if strcmp(handles.callerTag, 'selectDenoisingPopup')
        handles.asked_params = func([], [1 2 3]);
        handles.answered = true;
    else
        handles.asked_params = func([]);
        handles.answered = true;
    end
catch
    % Function didn't respond
    % No information can be gathered about the function's parameters
end
% Populate dialog if function answered
if handles.answered
    if ~isempty(handles.asked_params)
        
        % Check the data type of current selection and change interface if
        % different kind of input element is needed (editText for scalar
        % and pushButton for matrix).
        index_selected = get(handles.parameterPopup, 'Value');
        checkDataType(handles, index_selected);
        
        % --- Populate popUp
        if isfield(handles.asked_params, 'param')
            handles.numberOfParams = length(handles.asked_params.param);
            for i = 1:handles.numberOfParams
                itemToPopup(handles.parameterPopup, handles.asked_params.param{i});
            end
            % Initialize handles.params -struct with correct field names
            for i = 1:handles.numberOfParams
                handles.params.(handles.asked_params.param{i}) = [];
            end
        end
        
        % --- Populate editText
        
        % Case 1 mainGUI's current denoise/prepro/ortho function is different than what
        % is handled here --> Ask default parameter values from the
        % function itself
        if ~strcmp(func_name, func2str(handles.mainGUIhandles.params.(handles.func_abbr).h))
            % There must be 'param_value' field
            if isfield(handles.asked_params, 'param_value')
                handles.numberOfValues = length(handles.asked_params.param_value);
                % It must have at least one value
                if handles.numberOfValues > 0
                    set(handles.valueEditText, 'String', num2str(handles.asked_params.param_value{1}));
                end
            end
            % Fill handles.params -struct with values from asked parameters
            for i = 1:handles.numberOfParams
                if i > handles.numberOfValues
                    handles.params.(handles.asked_params.param{i}) = [];
                else
                    handles.params.(handles.asked_params.param{i}) = handles.asked_params.param_value{i};
                end
            end
            
            % Mark if there is a loaded matrix or struct parameter
            if handles.numberOfParams >= index_selected
                if ~isempty(handles.params.(handles.asked_params.param{index_selected}))
                    if isstruct(handles.params.(handles.asked_params.param{index_selected}))
                        set(handles.vectorText, 'String', 'Struct loaded.');
                    else
                        [M, N] = size(handles.params.(handles.asked_params.param{index_selected}));
                        set(handles.vectorText, 'String', ['Variable loaded (' num2str(M) 'x' num2str(N) ').']);
                    end
                else
                    set(handles.vectorText, 'String', 'Nothing loaded.');
                end
            end
        end
        
        % Case 2 mainGUI's current denoise/prepro/ortho function is the same as handled
        % here --> Ask current parameter values from mainGUI
        if strcmp(func_name, func2str(handles.mainGUIhandles.params.(handles.func_abbr).h))
            % There must be 'params' field
            if isfield(handles.mainGUIhandles.params.(handles.func_abbr), 'params')
                if ~isempty(handles.mainGUIhandles.params.(handles.func_abbr).params)
                    set(handles.valueEditText, 'String', num2str(handles.mainGUIhandles.params.(handles.func_abbr).params.(handles.asked_params.param{1})));
                end
            end
            
            % Fill handles.params -struct with values from mainGUI
            for i = 1:handles.numberOfParams
                handles.params.(handles.asked_params.param{i}) = handles.mainGUIhandles.params.(handles.func_abbr).params.(handles.asked_params.param{i});
                handles.numberOfValues = handles.numberOfValues + 1;
            end
        end
        
        % Mark if there is a loaded matrix or struct parameter
        if handles.numberOfParams >= index_selected
            if ~isempty(handles.params.(handles.asked_params.param{index_selected}))
                if isstruct(handles.params.(handles.asked_params.param{index_selected}))
                    set(handles.vectorText, 'String', 'Struct loaded.');
                else
                    [M, N] = size(handles.params.(handles.asked_params.param{index_selected}));
                    set(handles.vectorText, 'String', ['Variable loaded (' num2str(M) 'x' num2str(N) ').']);
                end
            else
                set(handles.vectorText, 'String', 'Nothing loaded.');
            end
        end
        
        % --- Populate infoText
        if isfield(handles.asked_params, 'param_desc')
            numberOfDescriptions = length(handles.asked_params.param_desc);
            if numberOfDescriptions > 0
                set(handles.infoText, 'String', handles.asked_params.param_desc{1});
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dss_gui_funcParams wait for user response (see UIRESUME)
% uiwait(handles.paramsGUI);


% --- Outputs from this function are returned to the command line.
function varargout = dss_gui_funcParams_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in parameterPopup.
function parameterPopup_Callback(hObject, eventdata, handles)
% hObject    handle to parameterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns parameterPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parameterPopup
function_strings = get(hObject, 'String');
index_selected = get(hObject, 'Value');

if handles.answered
    
    % Check the datatype of the parameter
    checkDataType(handles, index_selected);
    
    if ~isempty(handles.asked_params)
        % Populate editText with current value if there is one
        if strcmp(get(handles.valueEditText, 'Enable'), 'on')
            try
                set(handles.valueEditText, 'String', num2str(handles.params.(handles.asked_params.param{index_selected})));
            catch
                set(handles.valueEditText, 'String', '');
            end
        else
           % EditText is disabled, so we have a matrix/struct variable
           if handles.numberOfParams >= index_selected
                if ~isempty(handles.params.(handles.asked_params.param{index_selected}))
                    if isstruct(handles.params.(handles.asked_params.param{index_selected}))
                        set(handles.vectorText, 'String', 'Struct loaded.');
                    else
                        [M, N] = size(handles.params.(handles.asked_params.param{index_selected}));
                        set(handles.vectorText, 'String', ['Variable loaded (' num2str(M) 'x' num2str(N) ').']);
                    end
                else
                    set(handles.vectorText, 'String', 'Nothing loaded.');
                end
           end
        end
        % Populate infoText with description if there is one
        if isfield(handles.asked_params, 'param_desc')
            numberOfDescriptions = length(handles.asked_params.param_desc);
            if numberOfDescriptions >= index_selected
                if ~isempty(handles.asked_params.param_desc{index_selected})
                    set(handles.infoText, 'String', handles.asked_params.param_desc{index_selected});
                else
                    set(handles.infoText, 'String', 'No description available');
                end
            end
        end
    end
end


% --- Executes during object creation, after setting all properties.
function parameterPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function valueEditText_Callback(hObject, eventdata, handles)
% hObject    handle to valueEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valueEditText as text
%        str2double(get(hObject,'String')) returns contents of valueEditText as a double
index_selected = get(handles.parameterPopup, 'Value');
if length(handles.asked_params.param_type) >= index_selected
    if strcmp(handles.asked_params.param_type{index_selected}, 'function')
        if ~isempty(get(hObject, 'String'))
            handles.params.(handles.asked_params.param{index_selected}) = get(hObject, 'String');
        else
            handles.params = rmfield(handles.params, (handles.asked_params.param{index_selected}));
        end
    else
        if ~isempty(get(hObject, 'String'))
            handles.params.(handles.asked_params.param{index_selected}) = str2num(get(hObject, 'String'));
        else
            handles.params = rmfield(handles.params, (handles.asked_params.param{index_selected}));
        end
    end
else
    handles.params.(handles.asked_params.param{index_selected}) = str2num(get(hObject, 'String'));
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function valueEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valueEditText (see GCBO)
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
global LAST_DENOISING_FUNC;
global LAST_WHITENING_FUNC;
global LAST_ORTHO_FUNC;

% Only load parameters if function told its parameters
if handles.answered

    % Update parameters for appropriate function
    switch handles.callerTag
        case 'selectDenoisingPopup'
            if isfield(handles.mainGUIhandles.params.denf, 'params')
                handles.mainGUIhandles.params.denf = rmfield(handles.mainGUIhandles.params.denf, 'params');
            end
            % Disable beta function if previously one was selected
            if isfield(handles.mainGUIhandles.params, 'betaf')
                handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'betaf');
            end
            % If denoising function is tanh, select local beta as beta
            % function
            if strcmp(func2str(handles.func), 'denoise_tanh')
               handles.mainGUIhandles.params.betaf.h = @beta_tanh; 
            end
            handles.mainGUIhandles.params.denf.params = handles.params;
            handles.mainGUIhandles.params.denf.h = handles.func;
            LAST_DENOISING_FUNC = handles.valueOfPopup;
        case 'selectWhiteningPopup'
            if isfield(handles.mainGUIhandles.params.preprocf, 'params')
                handles.mainGUIhandles.params.preprocf = rmfield(handles.mainGUIhandles.params.preprocf, 'params');
            end
            handles.mainGUIhandles.params.preprocf.params = handles.params;
            handles.mainGUIhandles.params.preprocf.h = handles.func;
            LAST_WHITENING_FUNC = handles.valueOfPopup;
        case 'selectOrthoPopup'
            if isfield(handles.mainGUIhandles.params.orthof, 'params')
                handles.mainGUIhandles.params.orthof = rmfield(handles.mainGUIhandles.params.orthof, 'params');
            end
            handles.mainGUIhandles.params.orthof.params = handles.params;
            handles.mainGUIhandles.params.orthof.h = handles.func;
            LAST_ORTHO_FUNC = handles.valueOfPopup;
    end
    guidata(handles.mainGUI, handles.mainGUIhandles);
else
    switch handles.callerTag
        case 'selectDenoisingPopup'
            % Disable beta function if previously one was selected
            if isfield(handles.mainGUIhandles.params, 'betaf')
                handles.mainGUIhandles.params = rmfield(handles.mainGUIhandles.params, 'betaf');
            end
            handles.mainGUIhandles.params.denf.params = '';
            handles.mainGUIhandles.params.denf.h = handles.func;
            LAST_DENOISING_FUNC = handles.valueOfPopup;
        case 'selectWhiteningPopup'
            handles.mainGUIhandles.params.preprocf.params = '';
            handles.mainGUIhandles.params.preprocf.h = handles.func;
            LAST_WHITENING_FUNC = handles.valueOfPopup;
        case 'selectOrthoPopup'
            handles.mainGUIhandles.params.orthof.params = '';
            handles.mainGUIhandles.params.orthof.h = handles.func;
            LAST_ORTHO_FUNC = handles.valueOfPopup;
    end
    guidata(handles.mainGUI, handles.mainGUIhandles);
end

close(handles.paramsGUI);


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LAST_DENOISING_FUNC;
global LAST_WHITENING_FUNC;
global LAST_ORTHO_FUNC;
% Reselect previous selection if current is canceled
switch handles.callerTag
    case 'selectDenoisingPopup'
        set(handles.mainGUIhandles.selectDenoisingPopup, 'Value', LAST_DENOISING_FUNC);
    case 'selectWhiteningPopup'
        set(handles.mainGUIhandles.selectWhiteningPopup, 'Value', LAST_WHITENING_FUNC);
    case 'selectOrthoPopup'
        set(handles.mainGUIhandles.selectOrthoPopup, 'Value', LAST_ORTHO_FUNC);
end

close(handles.paramsGUI);


% --- Function for emptying a popUp menu
function cleanPopup(hObject)
% hObject   popUp menu to be emptied
% handles   handles to mainGUI objects
strings = {};
strings{1, 1} = 'Nothing available';
%strings{1, 2} = '---';
set(hObject, 'String', strings, 'Value', 1);


% --- Function for adding items to a popUp menu
%     Overwrites dummy item if there is one
function itemToPopup(hObject, item)
% hObject   popUp menu to receive item
% handles   handles to mainGUI objects
strings = get(hObject, 'String');
numberOfStrings = length(strings);
if strcmp(char(strings{1, 1}), 'Nothing available')
    strings{1, 1} = item;
    set(hObject, 'String', strings, 'Value', 1);
else
    strings{numberOfStrings+1, 1} = item;
    set(hObject, 'String', strings, 'Value', 1);
end


% --- Executes on button press in browseButton.
function browseButton_Callback(hObject, eventdata, handles)
% hObject    handle to browseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dss_gui_browse(hObject, handles, handles.paramsGUI, get(handles.parameterPopup, 'Value'), [],[]);


% --- Function for changing the type of input element (editText/Button)
%     based on what is the datatype of the current parameter
%     (scalar/matrix).
%     Used by checkDataType().
function changeInputType(handles)

% Check which element is visible and act accordingly
if strcmp(get(handles.browseButton, 'Enable'), 'off')
    % Make browseButton active, disable editText
    set(handles.browseButton, 'Enable', 'on', 'Visible', 'on');
    set(handles.valueEditText, 'Enable', 'off', 'Visible', 'off');
    set(handles.plotButton, 'Visible', 'on', 'Enable', 'on');
    set(handles.vectorText, 'Visible', 'on');
else
    % Make editText active, disable browseButton
    set(handles.browseButton, 'Enable', 'off', 'Visible', 'off');
    set(handles.valueEditText, 'Enable', 'on', 'Visible', 'on');
    set(handles.plotButton, 'Visible', 'off', 'Enable', 'on');
    set(handles.vectorText, 'Visible', 'off');
end


% --- Function for checking if input element should be changed.
%     Uses changeInputType() if change is needed.
function checkDataType(handles, index_selected)

% Check the datatype of the parameter
    if isfield(handles.asked_params, 'param_type')
        numberOfDataTypeFields = length(handles.asked_params.param_type);
        if numberOfDataTypeFields >= index_selected
            if ~isempty(handles.asked_params.param_type{index_selected})
                dataType = handles.asked_params.param_type{index_selected};
                switch dataType
                    case 'scalar'
                        set(handles.variableTypeText, 'String', 'scalar');
                        if strcmp(get(handles.browseButton, 'Enable'), 'on')
                            changeInputType(handles);
                        end
                    case 'vector'
                        set(handles.variableTypeText, 'String', 'vector');
                        if strcmp(get(handles.browseButton, 'Enable'), 'off')
                            changeInputType(handles);
                        end
                    case 'function'
                        set(handles.variableTypeText, 'String', 'function');
                        if strcmp(get(handles.browseButton, 'Enable'), 'on')
                            changeInputType(handles);
                        end
                    case 'structure'
                        set(handles.variableTypeText, 'String', 'structure');
                        if strcmp(get(handles.browseButton, 'Enable'), 'off')
                            changeInputType(handles);
                        end
                end 
            else
                setUnknownDataType(handles);
            end
        else
            setUnknownDataType(handles);
        end
    else
        setUnknownDataType(handles);
    end
    

% --- Function for setting datatype and input element when datatype is
%     unknown
function setUnknownDataType(handles)

% If datatype is not specified, use the editText
set(handles.variableTypeText, 'String', 'unknown');
if strcmp(get(handles.browseButton, 'Enable'), 'on')
    changeInputType(handles);
end


% --- Executes on button press in plotButton.
function plotButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

index_selected = get(handles.parameterPopup, 'Value');
if ~isempty(handles.params.(handles.asked_params.param{index_selected}))
    signal = handles.params.(handles.asked_params.param{index_selected});
    if isstruct(signal)
        return;
    end
    dim = size(signal, 1);
    [M, N] = size(signal);
    % Check if dim has a good value
    if dim < 21 | ( M > 1 & N > 1 )
        fig = figure; clf;
        if M > 1 & N > 1
            % We have a matrix
            set(fig, 'Name', 'Matrix parameter');
            imagesc(signal);
            colorbar;
        else
            set(fig, 'Name', 'Vector(s)');
            plotindex = 1;
            % Plot
            for sub_signal = 1:dim
                subplot(dim, 1, plotindex);
                plot(signal(sub_signal, :));
                yscale = max(-min(signal(sub_signal,:)), max(signal(sub_signal,:)))*1.3;
                axis([0 size(signal,2) -yscale yscale]);
                plotindex = plotindex + 1;
            end
        end
    else
        response = dss_gui_modalDlg('Title', ['Confirm plotting of ' num2str(dim) ' signals'], ...
            'String', 'Plotting that many signals may produce poor results. Are you sure?');
        switch response
            case 'Yes'
                set(fig, 'Name', 'Vector(s)');
                plotindex = 1;
                % Plot
                for sub_signal = 1:dim
                    subplot(dim, 1, plotindex);
                    plot(signal(sub_signal, :));
                    yscale = max(-min(signal(sub_signal,:)), max(signal(sub_signal,:)))*1.3;
                    axis([0 size(signal,2) -yscale yscale]);
                    plotindex = plotindex + 1;
                end
            case 'No'
                return;
        end
    end
end