function varargout = guireadneurone(varargin)
% GUIREADNEURONE Application M-file for guireadneurone.fig
%    FIG = GUIREADNEURONE launch guireadneurone GUI.
%    GUIREADNEURONE('callback_name', ...) invoke the named callback.
%
% ========================================================================
% NOTE:
% This file is part of the NeurOne data import plugin for EEGLAB.
% ========================================================================
% 
% Current version: 1.0.3.4 (2016-06-17)
% Author: Mega Electronics

% If no input arguments are used, the GUI is launched
if nargin == 0  
    
	neuroneimportfig = openfig(mfilename,'new');
    
	% Initialize a structure of handles to pass to callbacks.
	handles = guihandles(neuroneimportfig);
    
    % Load logos
    bgcolor=[0.656 0.758 1.0];
    megaLogo=imread('mega_gradient_edit.png','BackgroundColor',bgcolor);
    axes(handles.mega_logo);
    image(megaLogo) 
    axis off
    axis image
    
    neuroneLogo=imread('neurone_logo.png','BackgroundColor',bgcolor);
    axes(handles.neurone_logo);
    image(neuroneLogo)
    axis off
    axis image
    
    % Declare variables
    handles.chans='';
    handles.sessionPhaseNumber=1;
    handles.loadStatus=0;
    
    % Store the structure
	guidata(neuroneimportfig, handles)

	% Wait for callbacks. 'Ok' or 'Cancel' to continue.
	uiwait(neuroneimportfig);

	if nargout > 0
		varargout{1} = neuroneimportfig;
	end

elseif ischar(varargin{1})

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:});
		else
			feval(varargin{:});
		end
	catch
		disp(lasterr);
	end

end


% --- Executes on button press in cancel_button.
function varargout = cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.loadStatus=0;
guidata(hObject,handles);
uiresume(handles.guireadneurone_fig)


function varargout = channel_area_Callback(hObject, eventdata, handles)
% hObject    handle to channel_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel_area as text
%        str2double(get(hObject,'String')) returns contents of channel_area as a double

handles.chans=get(hObject,'string');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function varargout = channel_area_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varargout = session_area_Callback(hObject, eventdata, handles)
% hObject    handle to channel_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(h,'String') returns contents of channel_area as text
%        str2double(get(hObject,'String')) returns contents of channel_area as a double

handles.sessionPhaseNumber=str2num(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function varargout = session_area_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help_button.
function varargout = help_button_Callback(hObject, eventdata, handles)
% hObject    handle to help_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pophelp('pop_readneurone.m');


% --- Executes on button press in ok_button.
function varargout = ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
 handles.loadStatus=1;
 guidata(hObject,handles);
 uiresume(handles.guireadneurone_fig);
