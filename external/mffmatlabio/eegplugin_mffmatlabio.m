% eegplugin_mffmatlabio() - plugin for importing and exporting MFF files
%
% Usage:
%   >> eegplugin_mffmatlabio(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Main files
% ----------
% eegplugin_mffmatlabio.m - create menu in the EEGLAB interface
% pop_mffimport.m     - import MFF file into EEGLAB (graphic interface)
% pop_mffexport.m     - export EEGLAB structure to MFF file/folder (graphic interface)
%
% Support files
% -------------
% mff_import.m          - import MFF file/folder to EEGLAB structure
% mff_importcategories.m   - import 'categories.xml' file
% mff_importcoordinates.m  - import 'coordinates.xml' file
% mff_importepochs.m       - import 'epochs.xml' file
% mff_importevents.m       - import 'eventsxxxx.xml' file(s)
% mff_importinfo.m         - import 'info.xml' file
% mff_importinfon.m        - import 'info1.xml' file
% mff_importpnsset.m       - import PNS file
% mff_importsensorlayout.m - import 'sensorlayout.xml' file
% mff_importsignal.m       - import 'signal1.bin' file
% mff_importsubject.m      - import subject file
% mff_export.m             - export MFF file/folder from EEGLAB structure
% mff_createmff.m          - create empty MFF file/folder
% mff_exportcategories.m   - export 'categories.xml' file
% mff_exportcoordinates.m  - export 'coordinates.xml' file
% mff_exportepochs.m       - export 'epochs.xml' file
% mff_exportevents.m       - export 'events.xml' file
% mff_exportinfo.m         - export 'info.xml' file
% mff_exportinfon.m        - export 'info1.xml' file
% mff_exportpnsset.m       - export PNS file
% mff_exportsensorlayout.m - export 'sensorlayout.xml' file
% mff_exportsignal.m       - export 'signal1.bin' file
% mff_exportsubject.m      - export subject file
% mff_decodetime.m         - convert MFF time to Matlab time
% mff_encodetime.m         - convert Matlab time to MFF time
% mff_getobj.m             - convert Java object to Matlab
% mff_setobj.m             - convert Matlab object to Java
% pop_mffimport.m          - EEGLAB import GUI
% pop_mffexport.m          - EEGLAB export GUI
% eegplugin_mffmatlabio.m  - EEGLAB startup function
% eeg_compare.m            - Function to compare EEGLAB sturtures
% mff_fileio_read_data.m   - File-IO function to read data
% mff_fileio_read_header.m - File-IO function to read header
% mff_fileio_read_event.m  - File-IO function to read event
% mff_fileio_write.m       - File-IO function to write data and events

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function versionstr = eegplugin_mffmatlabio(fig, trystrs, catchstrs)

    %global EEG
    versionstr = '0.96';
    if nargin < 3
        error('eegplugin_mffmatlabio requires 3 arguments');
    end;
    
    % add amica folder to path
    % -----------------------
    if ~exist('mff_import')
        p = which('eegplugin_mffmatlabio');
        p = p(1:findstr(p,'eegplugin_mffmatlabio.m')-1);
        addpath(p);
    end

    % find tools menu
    % ---------------
    menui = findobj(fig, 'tag', 'import data'); 
    menue = findobj(fig, 'tag', 'export'); 
    
    % menu callback commands
    % ----------------------
    comload    = [  trystrs.no_check '[EEG, LASTCOM] = pop_mffimport;' catchstrs.store_and_hist ];
    comwrite   = [  trystrs.no_check 'LASTCOM = pop_mffexport(EEG);' catchstrs.store_and_hist ];
    
    % create menus
    % ------------
    submenu = uimenu( menui, 'Label', 'Import Philips .mff file', 'separator', 'on', 'CallBack', comload);
    submenu = uimenu( menue, 'Label', 'Export Philips .mff file', 'CallBack', comwrite);
