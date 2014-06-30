% icadefs() - function to read in a set of EEGLAB system-wide (i.e. lab-wide)
%             or working directory-wide constants and preferences. Change the 
%             way these are defined in the master icadefs.m file (usually
%             in dir eeglab/functions/sigprocfunc) or make a custom copy of 
%             the icadefs.m file in a project directory. Then, calling functions 
%             that call icadefs from an EEGLAB session in that working directory 
%             will read the local copy, which may set preferences different from 
%             the system-wide copy.
%
% Author: Arnaud Delorme, Scott Makeig, SCCN/INC/UCSD, La Jolla, 05-20-97 

% Copyright (C) 05-20-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% ----------------------------------------------------------------------
% ------ EEGLAB DEFINITION - YOU MAY CHANGE THE TEXT BELOW -------------
% ----------------------------------------------------------------------

EEGOPTION_PATH = ''; % if empty, the home folder of the current user is used
                     % Note that this may create problems under Windows
                     % when unicode characters are part of the user name
                     % In this case, enter the path name manually here.

YDIR  = 1;                  % positive potential up = 1; negative up = -1 
                            % for most ERP plots

HZDIR = 'up';               % ascending freqs = 'up'; descending = 'down' 
                            % (e.g., timef/newtimef frequency direction)

% font size
tmpComputer   = computer;
tmpScreenSize = get(0, 'ScreenSize');
retinaDisplay = false;
if tmpScreenSize(3) == 1440 && ( tmpScreenSize(3) == 878 || tmpScreenSize(3) == 900 )
    retinaDisplay = true;
end;
    
% retinaDisplay = false; % uncoment this line if not retina display
if retinaDisplay && strcmpi(tmpComputer(1:3), 'MAC')
    W_MAIN = findobj('tag', 'EEGLAB');
    if isempty(W_MAIN)
        disp('Mac OSX retina display detected. If this is not the case uncoment line 50 of icadefs.m');
    end;
    GUI_FONTSIZE  = 18; % graphic interface font size
    AXES_FONTSIZE = 18; % Axis labels and legend font size
    TEXT_FONTSIZE = 18; % Miscellaneous font sizes
else
    GUI_FONTSIZE  = 10; % graphic interface font size
    AXES_FONTSIZE = 10; % Axis labels and legend font size
    TEXT_FONTSIZE = 10; % Miscellaneous font sizes
end;
clear retinaDisplay tmpScreenSize tmpComputer;

% the eeg_options.m file also countains additional options

% ----------------------------------------------------------------------
% ------------------------ END OF DEFINITIONS --------------------------
% ----------------------------------------------------------------------

% INSERT location of ica executable (LINUX ONLY) for binica.m below
eeglab_p = fileparts(which('eeglab'));
ICABINARY = fullfile(eeglab_p, 'functions', 'resources', 'ica_linux'); 

try
    set(0,'defaultaxesfontsize',AXES_FONTSIZE);
    set(0,'defaulttextfontsize',TEXT_FONTSIZE);
    set(0,'DefaultUicontrolFontSize',GUI_FONTSIZE);
catch
    % most likely Octave here
end;

TUTORIAL_URL = 'http://sccn.ucsd.edu/wiki/EEGLAB'; % online version
DEFAULT_SRATE = 256.0175;      % default local sampling rate (rarely used)
DEFAULT_TIMLIM = [-1000 2000]; % default local epoch limits (ms)

% Set EEGLAB figure and GUI colors
% --------------------------------
lowscreendepth = 0;
if ~exist('OCTAVE_VERSION')
    if get(0, 'screendepth') <=8 % if mono or 8-bit color
	lowscreendepth = 1; 
    end;
end;
if lowscreendepth    
    %fprintf('icadefs(): Setting display parameters for mono or 8-bit color\n');
    BACKCOLOR           = [1 1 1];    % Background figure color 
    BACKEEGLABCOLOR     = [1 1 1];    % EEGLAB main window background
    GUIBUTTONCOLOR      = [1 1 1];    % Buttons colors in figures
    GUIPOPBUTTONCOLOR   = [1 1 1];    % Buttons colors in GUI windows
    GUIBACKCOLOR        = [1 1 1];    % GUI background color
    GUITEXTCOLOR        = [0 0 0];      % GUI foreground color for text    
    PLUGINMENUCOLOR     = [.5 0 .5];  % plugin menu color

else % if full color screen
    BACKCOLOR           = [.93 .96 1];    % EEGLAB Background figure color 
    BACKEEGLABCOLOR     = [.66 .76 1];    % EEGLAB main window background
    GUIBUTTONCOLOR      = BACKEEGLABCOLOR;% Buttons colors in figures
    GUIPOPBUTTONCOLOR   = BACKCOLOR;      % Buttons colors in GUI windows
    GUIBACKCOLOR        = BACKEEGLABCOLOR;% EEGLAB GUI background color <---------
    GUITEXTCOLOR        = [0 0 0.4];      % GUI foreground color for text
    PLUGINMENUCOLOR     = [.5 0 .5];      % plugin menu color
end;


% THE FOLLOWING PARAMETERS WILL BE DEPRECATED IN LATER VERSIONS
% -------------------------------------------------------------

SHRINKWARNING = 1;          % Warn user about the shrink factor in topoplot() (1/0)

MAXENVPLOTCHANS   = 264;  % maximum number of channels to plot in envproj.m
MAXPLOTDATACHANS  = 264;  % maximum number of channels to plot in dataplot.m
MAXPLOTDATAEPOCHS = 264;  % maximum number of epochs to plot in dataplot.m
MAXEEGPLOTCHANS   = 264;  % maximum number of channels to plot in eegplot.m
MAXTOPOPLOTCHANS  = 264;  % maximum number of channels to plot in topoplot.m

DEFAULT_ELOC  = 'chan.locs'; % default electrode location file for topoplot.m
DEFAULT_EPOCH = 10;       % default epoch width to plot in eegplot(s) (in sec)

SC  =  ['binica.sc'];           % Master .sc script file for binica.m
                                % MATLAB will use first such file found
                                % in its path of script directories.
                                % Copy to pwd to alter ICA defaults
