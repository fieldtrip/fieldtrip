function varargout = dss_gui_about(varargin)
% DSS_GUI_ABOUT M-file for dss_gui_about.fig
%      DSS_GUI_ABOUT, by itself, creates a new DSS_GUI_ABOUT or raises the existing
%      singleton*.
%
%      H = DSS_GUI_ABOUT returns the handle to a new DSS_GUI_ABOUT or the handle to
%      the existing singleton*.
%
%      DSS_GUI_ABOUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSS_GUI_ABOUT.M with the given input arguments.
%
%      DSS_GUI_ABOUT('Property','Value',...) creates a new DSS_GUI_ABOUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dss_gui_about_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dss_gui_about_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Includes functionality for the about window.
%      This file is used by other DSS files and should not be used alone.
%
% See also: GUIDE, GUIDATA, DSS_GUI_BROWSE, DSS_GUI_ADVOPTIONS,
% DSS_GUI_FUNCPARAMS, DSS_GUI_INSERTFUNC, DSS_GUI_MODALDLG, DSS_GUI_SAVE,
% DSS_GUI_ABOUT, DSS_GUI_HELP

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss.

% Edit the above text to modify the response to help dss_gui_about

% Last Modified by GUIDE v2.5 28-Feb-2005 10:54:16

% $Id: dss_gui_about.m,v 1.7 2005/04/20 10:19:24 kosti Exp $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dss_gui_about_OpeningFcn, ...
                   'gui_OutputFcn',  @dss_gui_about_OutputFcn, ...
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


% --- Executes just before dss_gui_about is made visible.
function dss_gui_about_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dss_gui_about (see VARARGIN)

% Choose default command line output for dss_gui_about
handles.output = hObject;

% Name dialog
set(hObject, 'Name', 'About', 'WindowStyle', 'modal');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dss_gui_about wait for user response (see UIRESUME)
% uiwait(handles.aboutGUI);


% --- Outputs from this function are returned to the command line.
function varargout = dss_gui_about_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.aboutGUI);