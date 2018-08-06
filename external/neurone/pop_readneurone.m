% pop_readneurone() - A function for reading NeurOne data.
%
% Usage: >> [EEGOUT,command] = pop_readneurone() % Launches GUI for user input
%       OR
%        >> [EEGOUT,command] = pop_readneurone(dataPath, sessionPhaseNumber, chans)
%
% Please note that if this function is called without using the menu item
% in EEGLAB, it will not be automatically stored in its memory.
% ========================================================================
% Inputs:
%        dataPath           - A string containing the full path for .ses
%                             -file to be read. This input argument is
%                             required.
%        sessionPhaseNumber - A number indicating the session
%                             phase to be observed. This input argument is
%                             optional. If the session phase number is left
%                             empty, default value 1 will be used.
%        chans              - Another optional input argument to indicate
%                             the channels to be imported. Acceptable
%                             values are strings where the channels are 
%                             in a vector form without brackets e.g. 
%                             '1,3,5' or '1:4' (when the GUI is executed, 
%                             the quotation marks should not be used). 
%                             Note that only acceptable values
%                             are those between 1 and the total number of
%                             channels. It should also be noted that the
%                             number of the channel may not be the same as
%                             its input number since when an input
%                             number is missing, the channel numbering will
%                             continue with next available input number.
%                             Here is an example of the channel numbering
%                             (note the missing input numbers 3 and 4):
%
%              Input number               Channel number
%             ------------               --------------
%                   1          ----->           1
%                   2          ----->           2
%                   5          ----->           3
%                   6          ----->           4
%
%
% Outputs:
%        EEGOUT             - A struct with the same fields as a standard 
%                             EEG data structure under EEGLAB 
%                             (see eeg_checkset.m).
%        command            - A string containing all the input parameters
%                             and the command used to execute this
%                             function.
% ========================================================================
% GUI help:
%
% When this function is used without any input arguments, it will open
% several windows to ask for required information as described above. The
% first pop-up window asks the user to specify the .ses -file to be loaded.
% If 'Cancel' -button is pressed at this point, the importing process will
% terminate. The next window is a GUI asking for information about the
% channels to be loaded and the number of the desired session phase to be
% observed. When 'Ok' -button is pressed, the data is sent to readneurone.m
% for further processing. Again, if 'Cancel' -button is pressed, the whole
% importing process will be terminated. 
% 
% IMPORTANT NOTE: For a computer with high memory capacity it is much more
%                 faster to import data by leaving the channel field empty
%                 and removing unnecessary channels afterwards. However,
%                 in case of limited memory it is recommended to first
%                 specify the channels to be imported.
%
% ========================================================================
% NOTE:
% This file is part of the NeurOne data import plugin for EEGLAB.
% ========================================================================
% 
% Current version: 1.0.3.4 (2016-06-17)
% Author: Mega Electronics
%
% see also: eegplugin_neurone.m or readneurone.m

function [EEG,command] = pop_readneurone(dataPath,sessionPhaseNumber,chans)
%% Initialize empty structure for the data

EEG = {};              % Create empty space for data
command = '';

%% Check the number of input arguments
nargchk(0,1,nargin);

%% Use the GUI

if nargin==0
    
    if isfield(evalin('base','EEG'),'history')
        history = evalin('base','EEG.history');
        if  ~(isempty(strfind(history,'dataPath')))
            dataPathIndex = strfind(history,'dataPath=''');
            idx = strfind(history,'''');
            idx = idx(idx>dataPathIndex);
            dataPath = history(idx(1)+1:idx(2)-1);
            dataPath = fileparts(dataPath);
            [dataPath folder] = fileparts(dataPath);
        else
            dataPath = cd;
        end
    else
        dataPath = cd;
    end
    
    [fname dataPath] = uigetfile({'*.ses', ...
        'NeurOne session file (*.ses)'}, ...
        'Load NeurOne session file',dataPath);
    
    % If cancel is pressed, the import process will terminate
    if fname==0
        return
    end
    
    % If ok, manipulate the filename
    fullpath = [dataPath fname];
    a = dir(dataPath);
    for k = 3:numel(a)
        if a(k).isdir==1
            if ~(isempty(regexpi(fullpath,a(k).name)))
                dataPath = [dataPath a(k).name filesep];
            end
        end
    end
    
    fprintf('\n----------IMPORTING DATA----------\n')
    fprintf('File to Import: %s\n',fullpath)
      
    guifig=guireadneurone;

    try
       SETTINGS=guidata(guifig);
    catch
       return;
    end
    
    switch SETTINGS.loadStatus
    % When cancel button is pressed, the GUI terminates
    case 0
        disp('NeurOne import cancelled')
        delete(guifig);
        return
    % When ok button is pressed, proceed to read the data
    case 1
        delete(guifig);   

        SETTINGS.dataPath=dataPath;

        if  isempty(SETTINGS.chans)
            fprintf('Channels to be read: all\n');
        else
            fprintf('Channels to be read: %s\n',SETTINGS.chans);
        end

        EEG = readneurone(SETTINGS.dataPath, ...
            SETTINGS.sessionPhaseNumber,SETTINGS.chans);
        EEG = eeg_checkset(EEG); % EEGLAB function to check consistency
        EEG = eeg_checkset(EEG,'eventconsistency'); % check event structure
        EEG = eeg_checkset(EEG,'makeur');
    end

% ====================================================================
% Read data using direct command
    
else
    if nargin>=1
         fullpath=dataPath;
         fprintf('\n----------IMPORTING DATA----------\n')
         fprintf('File to Import: %s\n',fullpath)
    
      % Manipulate the given path
      [dataPath fname] = fileparts(fullpath);
      a = dir(dataPath);
      for k = 3:numel(a)
          if a(k).isdir==1
             if ~(isempty(regexpi(fullpath,a(k).name)))
                dataPath = [dataPath filesep a(k).name filesep];
             end
          end
      end  
         SETTINGS=struct();
         SETTINGS.dataPath=dataPath;
         
         if ~(exist('sessionPhaseNumber'))
            SETTINGS.sessionPhaseNumber = 1;
         else
            SETTINGS.sessionPhaseNumber=sessionPhaseNumber;
         end
         
         if ~(exist('chans'))
            SETTINGS.chans = '';
         else
            SETTINGS.chans=chans;
         end
        
    end
    
    EEG = readneurone(SETTINGS.dataPath,SETTINGS.sessionPhaseNumber, ...
        SETTINGS.chans);
    EEG = eeg_checkset(EEG);                    % check EEG data structure
    EEG = eeg_checkset(EEG,'eventconsistency'); % check event structure
    EEG = eeg_checkset(EEG,'makeur');

end

comline = sprintf('dataPath=''%s'';',SETTINGS.dataPath);
comline = sprintf('%s sessionPhaseNumber=''%d'';',...
    comline,SETTINGS.sessionPhaseNumber);
comline = sprintf('%s chans=''%s'';',...
    comline,SETTINGS.chans);
command = sprintf('%s [EEG,COM]=pop_readneurone(dataPath,sessionPhaseNumber,chans);',comline);  


return