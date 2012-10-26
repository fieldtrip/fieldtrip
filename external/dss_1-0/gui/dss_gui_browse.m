function varargout = dss_gui_browse(varargin)
% DSS_GUI_BROWSE M-file for dss_gui_browse.fig
%      DSS_GUI_BROWSE, by itself, creates a new DSS_GUI_BROWSE or raises the existing
%      singleton*.
%
%      H = DSS_GUI_BROWSE returns the handle to a new DSS_GUI_BROWSE or the handle to
%      the existing singleton*.
%
%      DSS_GUI_BROWSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSS_GUI_BROWSE.M with the given input arguments.
%
%      DSS_GUI_BROWSE('Property','Value',...) creates a new DSS_GUI_BROWSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dss_gui_browse_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dss_gui_browse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Includes functionality of the browse window.
%      This file is used by other DSS files and should not be used alone.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% Edit the above text to modify the response to help dss_gui_browse

% Last Modified by GUIDE v2.5 16-Feb-2005 11:46:03

% $Id: dss_gui_browse.m,v 1.8 2005/04/20 10:19:24 kosti Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dss_gui_browse_OpeningFcn, ...
                   'gui_OutputFcn',  @dss_gui_browse_OutputFcn, ...
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


% --- Executes just before dss_gui_browse is made visible.
function dss_gui_browse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dss_gui_browse (see VARARGIN)
global DSS_DIRECTORY;

% Choose default command line output for dss_gui_browse
handles.output = hObject;

% Handles from main GUI
handles.funcParamIndex = varargin{4};
handles.mainGUI = varargin{3};
handles.mainGUIhandles = varargin{2};
handles.caller = varargin{1};
handles.stateDefaultsFcn = varargin{5};
handles.defaultsFcn = varargin{6};

% Find out who called this dialog
guidata(hObject, handles);
handles.callerString = get(handles.caller, 'String');

% Make dialog modal
set(hObject, 'WindowStyle','modal');

% Update handles structure
guidata(hObject, handles);

% Showing workspace or working directory based on who called
if strcmp(handles.callerString, 'Browse files')
    set(hObject, 'Name', 'Browse files');
    set(handles.text1, 'String', 'Files');
    
    handles.starting_dir = pwd;
     
    % Populate the listbox with files
    handles = refreshFileList(pwd, handles);
    showDimensions(handles);
elseif strcmp(handles.callerString, 'Browse')
    set(hObject, 'Name', 'Browse workspace');
    % Refresh the listbox (workspace)
    refreshVariables(handles);
    showDimensions(handles);
elseif strcmp(handles.callerString, 'Browse demos')
    set(hObject, 'Name', 'Browse demo files');
    set(handles.text1, 'String', 'Demo data');
    
    try
        handles.starting_dir = pwd;
        cd (DSS_DIRECTORY)
        cd demo
    catch
        errordlg('DSS demo-directory not found.', 'Directory Error', 'modal');
        cd (handles.starting_dir)
        close(hObject);
    end
    handles = refreshFileList(pwd, handles);
    showDimensions(handles);
elseif strcmp(handles.callerString, 'Load params')
    set(hObject, 'Name', 'Browse for function params struct');
    % Refresh the listbox (workspace)
    refreshVariables(handles);
    showDimensions(handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dss_gui_browse wait for user response (see UIRESUME)
% uiwait(handles.loadGUI);


% --- Outputs from this function are returned to the command line.
function varargout = dss_gui_browse_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.caller, 'Tag'), 'browseButton')
    loadFuncParam(handles, handles.funcParamIndex);
else
    loadData(handles);
end


% --- Executes on button press in plotButton.
function plotButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.callerString, 'Browse')
    signal_string = get_var_names(handles);
    signal = evalin('base', signal_string);
    if isempty(signal_string)
        signal = [];
    end
    if ~isempty(signal) & isnumeric(signal)
        dim = size(signal, 1);
        [M, N] = size(signal);
        % Check if dim has a good value
        if dim < 21 | ( strcmp(get(handles.caller, 'Tag'), 'browseButton') & M > 1 & N > 1 )
            fig = figure; clf;
            if strcmp(get(handles.caller, 'Tag'), 'browseButton') & M > 1 & N > 1
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
                    fig = figure; clf;
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
    else
        errordlg('Data must be non empty and numeric.', 'Input Error', 'modal');
    end
    
elseif strcmp(handles.callerString, 'Browse files') | strcmp(handles.callerString, 'Browse demos')
    % Browsing files
    signal_string = get_var_names(handles);
    [path, name, ext, ver] = fileparts(signal_string);
    % File must be .mat file
    if ~strcmp(ext, '.mat')
        errordlg('File must be .mat file', 'Input error', 'modal');
        return;
    elseif exist(signal_string) ~= 2
        errordlg('File does not exist.', 'Input error', 'modal');
        return;
    else
        % Load file
        tempvar = load(signal_string);
        varnames = fieldnames(tempvar);
        signal = tempvar.(varnames{1});
    end
    if ~isempty(signal) & isnumeric(signal)
        dim = size(signal, 1);
        % Check if dim has a good value
        if dim < 21
            fig = figure; clf;
            set(fig, 'Name', 'Signal');
            plotindex = 1;
            % Plot
            for sub_signal = 1:dim
                subplot(dim, 1, plotindex);
                plot(signal(sub_signal, :));
                yscale = max(-min(signal(sub_signal,:)), max(signal(sub_signal,:)))*1.3;
                axis([0 size(signal,2) -yscale yscale]);
                plotindex = plotindex + 1;
            end
        else
            response = dss_gui_modalDlg('Title', ['Confirm plotting of ' num2str(dim) ' signals'], ...
                'String', 'Plotting that many signals may produce poor results. Are you sure?');
            switch response
                case 'Yes'
                    fig = figure; clf;
                    set(fig, 'Name', 'Signal');
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
    else
        errordlg('Data must be non empty and numeric.', 'Input Error', 'modal');
    end
end


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.callerString, 'Browse demos') | strcmp(handles.callerString, 'Browse files')
    if isfield(handles, 'starting_dir')
        cd (handles.starting_dir)
    end
end
close(handles.loadGUI);


% --- Executes on button press in refreshButton.
function refreshButton_Callback(h, eventdata, handles, varargin)
% hObject    handle to refreshButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.callerString, 'Browse files') | strcmp(handles.callerString, 'Browse demos')
    refreshFileList(pwd, handles);
elseif strcmp(handles.callerString, 'Browse') | strcmp(handles.callerString, 'Load params')
    refreshVariables(handles);
end
showDimensions(handles);


function refreshVariables(handles)
% hObject    handle to update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% Updates the listbox to match the current workspace

vars = evalin('base','who');
set(handles.listbox1, 'String', vars);
selected_value = get(handles.listbox1,'Value');

% Fix for bad behaviour if user deselects everything in the list
if length(vars) > 0 & isempty(selected_value)
        set(handles.listbox1, 'Value', 1);
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% If double clicked and we are browsing files
if strcmp(get(handles.loadGUI, 'SelectionType'), 'open') & (strcmp(handles.callerString, 'Browse files') | strcmp(handles.callerString, 'Browse demos'))
    index_selected = get(handles.listbox1, 'Value');
    file_list = get(handles.listbox1, 'String');
    filename = file_list{index_selected};
    if handles.is_dir(handles.sorted_index(index_selected))
        cd (filename);
        refreshFileList(pwd, handles);
    else
        [path, name, ext, ver] = fileparts(filename);
        switch ext
            case '.mat'
                loadData(handles);
            otherwise
                errordlg('File must be .mat file.', 'File Type Error', 'modal');
        end
    end
elseif strcmp(get(handles.loadGUI, 'SelectionType'), 'open') & (strcmp(handles.callerString, 'Browse') | strcmp(handles.callerString, 'Load params'))
    if strcmp(get(handles.caller, 'Tag'), 'browseButton')
        loadFuncParam(handles, handles.funcParamIndex);
    else
        loadData(handles);
    end
elseif strcmp(get(handles.loadGUI, 'SelectionType'), 'normal')
    % Normal selection (one click)
    showDimensions(handles);
end


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% Function for getting the selected variable's name
function varname = get_var_names(handles)

list_entries = get(handles.listbox1, 'String');
index_selected = get(handles.listbox1, 'Value');
% Exactly one value must be selected
if length(index_selected) ~= 1
    if length(list_entries) > 0
        set(handles.listbox1, 'Value', 1);
    end
    varname = '';
else
    if ~isempty(list_entries)
        varname = list_entries{index_selected};
    else
        varname = '';
    end
end


% Function for listing files in current directory
function handles = refreshFileList(dir_path, handles)

cd (dir_path)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
temp_names = {};
k = 1;
% Sorting the names again, mat files to top
for i = 1:length(sorted_names)
    % Path to lower dir
    [path, name, ext, ver] = fileparts(sorted_names{i});
    if strcmp(ext, '.')
        temp_names{k} = sorted_names{i};
        sorted_index(k) = i;
        k = k + 1;
    end
end
for i = 1:length(sorted_names)
    % Mat files
    [path, name, ext, ver] = fileparts(sorted_names{i});
    if strcmp(ext, '.mat')
        temp_names{k} = sorted_names{i};
        sorted_index(k) = i;
        k = k + 1;
    end
end
for i = 1:length(sorted_names)
    % Directories
    [path, name, ext, ver] = fileparts(sorted_names{i});
    if strcmp(ext, '')
        temp_names{k} = sorted_names{i};
        sorted_index(k) = i;
        k = k + 1;
    end
end
for i = 1:length(sorted_names)
    % Everything else
    [path, name, ext, ver] = fileparts(sorted_names{i});
    if ~strcmp(ext, '.mat') & ~strcmp(ext, '.') & ~strcmp(ext, '')
        temp_names{k} = sorted_names{i};
        sorted_index(k) = i;
        k = k + 1;
    end
end

handles.file_names = temp_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = [sorted_index];
guidata(handles.loadGUI, handles);
set(handles.listbox1, 'String', handles.file_names, 'Value', 1);


% --- Function for loading data
function loadData(handles)

global INPUTSIGNAL;
global INPUTSIGNAL_STRING;

% Loading variable from workspace
if strcmp(handles.callerString, 'Browse')
    varname = get_var_names(handles);
    if ~isempty(varname)
        tempvar = evalin('base', varname);
        % Variable must be numeric
        if ~isnumeric(tempvar)
            if isstruct(tempvar)
                % Load state struct
                handles.mainGUIhandles.state_given = true;
                handles.mainGUIhandles.state = tempvar;
                INPUTSIGNAL = [];
                INPUTSIGNAL_STRING = varname;
                [signals, samples] = size(tempvar.X);
                refreshMainGUI(handles, signals, samples);
                handles.mainGUIhandles = feval(handles.stateDefaultsFcn, handles.mainGUIhandles, true, true);
                guidata(handles.mainGUI, handles.mainGUIhandles);
                close(handles.loadGUI);
            else
                errordlg('Variable must be numeric or struct.', 'Input error', 'modal');
                return;
            end
        elseif isempty(tempvar)
            errordlg('Variable is empty.', 'Input error', 'modal');
            return;
        else
            % Load variable
            INPUTSIGNAL = tempvar;
            INPUTSIGNAL_STRING = varname;
            
            [signals, samples] = size(INPUTSIGNAL);
            refreshMainGUI(handles, signals, samples);

            if handles.mainGUIhandles.state_given
                handles.mainGUIhandles.state_given = false;
                handles.mainGUIhandles = feval(handles.defaultsFcn, handles.mainGUIhandles);
            end
            guidata(handles.mainGUI, handles.mainGUIhandles);
            close(handles.loadGUI);
        end
    else
        errordlg('Nothing selected.', 'Input error', 'modal');
    end
% Loading a file    
elseif strcmp(handles.callerString, 'Browse files') | strcmp(handles.callerString, 'Browse demos')
    varname = get_var_names(handles);
    if ~isempty(varname)
        [path, name, ext, ver] = fileparts(varname);
        % File must be .mat file
        if ~strcmp(ext, '.mat')
            errordlg('File must be .mat file', 'Input error', 'modal');
            return;
        elseif exist(varname) ~= 2
            errordlg('File does not exist.', 'Input error', 'modal');
            return;
        else
            % Load file
            tempvar = load(varname);
            varnames = fieldnames(tempvar);
            var_count = length(varnames);
            
            if isstruct(tempvar.(varnames{1}))
                % Load state
                handles.mainGUIhandles.state_given = true;
                handles.mainGUIhandles.state = tempvar.(varnames{1});
                INPUTSIGNAL = [];
                INPUTSIGNAL_STRING = varname;
                [signals, samples] = size(tempvar.(varnames{1}).X);
                refreshMainGUI(handles, signals, samples);
                handles.mainGUIhandles = feval(handles.stateDefaultsFcn, handles.mainGUIhandles, true, true);
                guidata(handles.mainGUI, handles.mainGUIhandles);
                close(handles.loadGUI);
            else
                % Load signal
                INPUTSIGNAL = tempvar.(varnames{1});
                INPUTSIGNAL_STRING = varname;

                % TODO: loading of multiple variables in single file
                %{
                for i = 1:var_count
                    temp = tempvar.(varnames{i});
                    INPUTSIGNAL(i) = temp(:,:);
                end
                %}

                [signals, samples] = size(INPUTSIGNAL);
                refreshMainGUI(handles, signals, samples);

                if isfield(handles, 'starting_dir')
                    cd (handles.starting_dir)
                end
                if handles.mainGUIhandles.state_given
                    handles.mainGUIhandles.state_given = false;
                    handles.mainGUIhandles = feval(handles.defaultsFcn, handles.mainGUIhandles);
                end
                guidata(handles.mainGUI, handles.mainGUIhandles);
                close(handles.loadGUI);
            end
        end
    else
        errordlg('Nothing selected.', 'Input error', 'modal');
    end
end


% --- Function for loading a matrix as a function's parameter
function loadFuncParam(handles, index)

varname = get_var_names(handles);
if ~isempty(varname)
    tempvar = evalin('base', varname);
    if isstruct(tempvar)
        handles.mainGUIhandles.params.(handles.mainGUIhandles.asked_params.param{index}) = tempvar;
        handles.mainGUIhandles.varnames{index} = varname;
        guidata(handles.mainGUI, handles.mainGUIhandles);
        % Refresh paramsGUI
        set(handles.mainGUIhandles.vectorText, 'Visible', 'on', 'String', 'Struct loaded.');
        close(handles.loadGUI);
    else
        % Variable must be numeric
        if ~isnumeric(tempvar)
            errordlg('Variable must be numeric.', 'Input error', 'modal');
            return;
        elseif isempty(tempvar)
            errordlg('Variable is empty.', 'Input error', 'modal');
            return;
        else
            % Load variable
            [M, N] = size(tempvar);
            handles.mainGUIhandles.params.(handles.mainGUIhandles.asked_params.param{index}) = tempvar;
            handles.mainGUIhandles.varnames{index} = varname;
            guidata(handles.mainGUI, handles.mainGUIhandles);
            % Refresh paramsGUI
            set(handles.mainGUIhandles.vectorText, 'Visible', 'on', 'String', ['Variable loaded (' num2str(M) 'x' num2str(N) ').']);
            close(handles.loadGUI);
        end
    end
else
    errordlg('Nothing selected.', 'Input error', 'modal');
end


% --- Function for showing selected item's dimensions below the listbox
function showDimensions(handles)  
% Show the selected variables dimensions
var_name = get_var_names(handles);
if ~isempty(var_name)
    if strcmp(handles.callerString, 'Browse') | strcmp(handles.callerString, 'Load params')
        % We are browsing workspace
        var = evalin('base', var_name);
        if isnumeric(var)
            [M, N] = size(var);
            set(handles.infoText, 'String', ['Matrix: ' num2str(M) 'x' num2str(N)]);
        elseif isstruct(var)
            set(handles.infoText, 'String', 'Struct');
        else
            set(handles.infoText, 'String', 'Not numeric.');
        end
    else
        % We are browsing files
        [path, name, ext, ver] = fileparts(var_name);
        % File must be .mat file
        if ~strcmp(ext, '.mat') & ~strcmp(ext, '') & ~strcmp(ext, '.')
            set(handles.infoText, 'String', 'Unsupported file');
            return;
        elseif strcmp(ext, '') | strcmp(ext, '.')
            set(handles.infoText, 'String', 'Directory');
            return;
        end
        tempvar = load(var_name);
        varnames = fieldnames(tempvar);
        var = tempvar.(varnames{1});
        if isnumeric(var)
           [M, N] = size(var);
           set(handles.infoText, 'String', ['Matrix: ' num2str(M) 'x' num2str(N)]);
        elseif isstruct(var)
            set(handles.infoText, 'String', 'Struct');
        else
            set(handles.infoText, 'String', 'Not numeric.');
        end
    end
end


% --- Function for making some visual settings in main GUI, when new data
%     is loaded
function refreshMainGUI(handles, signals, samples)
% signals   number of signals in the new data
% samples   number of samples in each signal in the new data
global INPUTSIGNAL_STRING;

% Refresh main GUI
% Enable some elements
set(handles.mainGUIhandles.plotButton, 'Enable', 'on');
set(handles.mainGUIhandles.transposeButton, 'Enable', 'on');
set(handles.mainGUIhandles.startDssButton, 'Enable', 'on');
set(handles.mainGUIhandles.systemResetButton, 'Enable', 'on');
set(handles.mainGUIhandles.numberOfSignalsText, 'String', num2str(signals));
set(handles.mainGUIhandles.numberOfSamplesText, 'String', num2str(samples));
set(handles.mainGUIhandles.wdimInputBox, 'String', num2str(signals));
set(handles.mainGUIhandles.sdimEditText, 'String', num2str(signals));
if handles.mainGUIhandles.state_given
    set(handles.mainGUIhandles.informText, 'String', ['State ' INPUTSIGNAL_STRING ' loaded.']);
else
    set(handles.mainGUIhandles.informText, 'String', ['Signal ' INPUTSIGNAL_STRING ' loaded.']);
end