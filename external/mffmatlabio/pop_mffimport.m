% pop_mffimport - import MFF file to EEGLAB structure. Note that events
%                 in MFF files are very rich and that you need to specify
%                 which field may contain the information for EEGLAB to
%                 extract events. This is the purpose of the second input
%                 "typefield" to that function.
%
% Usage:
%   EEG = pop_mffimport; % pop up menu to select file
%   EEG = pop_mffimport(mffFile, typefield); % import file
%
% Input:
%  mffFile - filename/foldername for the MFF file (MFF file/folder must
%            already exist)
%  typefield - [string or cell] MFF field(s) to use for the type field. The
%            default is the code field of the MFF file. If several fields
%            are provided, the fields are concatenated to create the EEGLAB
%            type field.
%  savedat - [0|1] automatically save imported dataset
%
% Output:
%  EEG     - EEGLAB structure

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

function [EEG, com] = pop_mffimport(fileName, typefield, saveData)

com = '';
saveData = 0;
matVer = ver('MATLAB');
if datenum(matVer.Date) < 735595
    error('This version of Matlab is too old. Use version 2014a or later');
end

EEG = [];

if ~exist('convertlocs.m', 'file')
    error('This function requires to install the EEGLAB software (and start it)');
end

if nargin < 1
    
    % pop up window
    % -------------
    if ismac
        [fileName, filePath] = uigetfile('*', 'Select an EGI .mff file(s)', 'MultiSelect', 'on');
        if isnumeric(fileName) && fileName(1) == 0, return; end
        fileName = fullfile(filePath, fileName);
    else
        fileName = uigetdir2('*', 'Select an EGI .mff folder(s)');
        if isempty(fileName), return; end
        if length(fileName) == 1, fileName = fileName{1}; end
    end
    
    if iscell(fileName)
        buttonName = questdlg2([ 'Do you want to automatically save imported datasets?' 10 ...
            '(the name will remain the same as the original dataset' 10 ...
            'and the .set extension will be used)' ], 'pop_importmff() - import MFF dataset(s)', 'Cancel', 'No thanks', 'Yes Save', 'Yes Save');
        switch buttonName
            case 'Cancel', return;
            case 'No thanks', saveData = 0;
            otherwise saveData = 1;
        end
    end
end

if ~iscell(fileName), fileName = { fileName }; end
EEGTMP = [];
for iFile = 1:length(fileName)
    EEGTMP = mff_import(fileName{iFile});
    
    if nargin < 2 && iFile == 1
        if isempty(EEGTMP.event)
            disp('No event in data file');
        else
            % popup window parameters
            % -----------------------
            eventfields = setdiff(fieldnames(EEGTMP.event), { 'type' 'latency' 'duration' 'urevent' 'epoch' 'begintime' });
            poscode     = strmatch('code', lower(eventfields), 'exact');
            if isempty(poscode), poscode=1; end
            uilist   = { { 'style' 'text' 'String' ['Event type field:' 10 '(you may select multiple)'] } ...
                { 'style' 'listbox' 'string' eventfields 'value' poscode 'min' 0 'max' 2} };
            geom = { 1 1 };
            result = inputgui( 'geometry', geom, 'uilist', uilist, 'geomvert', [1 3], 'helpcom', 'pophelp(''pop_mffimport'')', 'title', 'Choose event type field -- pop_mffimport()');
            
            if isempty(result), return; end
            typefield = eventfields(result{1});
        end
    end
    
    % Use different event fields to populate the EEGLAB type field
    % ------------------------------------------------------------
    if exist('typefield', 'var')
        if ~isempty(typefield) || ~(ischar(typefield) && strcmpi(typefield, 'code'))
            if ischar(typefield), typefield = { typefield }; end
            
            % get data for each MFF event field
            dataField = cell(length(EEGTMP.event), length(typefield));
            for iField = 1:length(typefield)
                dataField(:,iField) = { EEGTMP.event.(typefield{iField}) }';
            end
            
            % copy the data into EEGLAB event structure
            strField = cell(1,length(typefield));
            strField(:) = {'_'}; % add space between events
            for iEvent = 1:length(EEGTMP.event)
                if ~isequal(lower(EEGTMP.event(iEvent).type), 'boundary') % not a boundary
                    if ~all(cellfun(@isempty, dataField(iEvent,:))) % not a trial type
                        tmpType = [ dataField(iEvent,:); strField ];
                        tmpType = [ tmpType{:} ];
                        tmpType(end) = [];
                        EEGTMP.event(iEvent).type = tmpType;
                    end
                end
            end
        end
    end
    EEGTMP = eeg_checkset(EEGTMP,'eventconsistency');
    if saveData
        EEGTMP = pop_saveset(EEGTMP, [ fileName{iFile}(1:end-4) '.set' ]);
    end
    EEG = eeg_store(EEG, EEGTMP);
end

% com = sprintf('EEG = pop_mffimport(''%s'', %s);', fileName, vararg2str({typefield}));
com = sprintf('EEG = pop_mffimport(%s', vararg2str({fileName}));
if exist('typefield', 'var'), com = sprintf([com  ',%s'],vararg2str({typefield})); end
if saveData, com = [ com ',1' ]; end
com = [com ');'];

% function below downlaoded from https://www.mathworks.com/matlabcentral/fileexchange/32555-uigetfile_n_dir-select-multiple-files-and-directories
% Copyright (c) 2011, Peugas
% All rights reserved.
function [pathname] = uigetdir2(start_path, dialog_title)
% Pick multiple directories and/or files

import javax.swing.JFileChooser;

if nargin == 0 || isempty(start_path)
    start_path = pwd;
elseif numel(start_path) == 1
    if start_path == 0 % Allow a null argument.
        start_path = pwd;
    end
end

jchooser = javaObjectEDT('javax.swing.JFileChooser', start_path);

jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
if nargin > 1
    jchooser.setDialogTitle(dialog_title);
end

jchooser.setMultiSelectionEnabled(true);

status = jchooser.showOpenDialog([]);

if status == JFileChooser.APPROVE_OPTION
    jFile = jchooser.getSelectedFiles();
	pathname{size(jFile, 1)}=[];
    for i=1:size(jFile, 1)
		pathname{i} = char(jFile(i).getAbsolutePath);
	end
	
elseif status == JFileChooser.CANCEL_OPTION
    pathname = [];
else
    error('Error occured while picking file.');
end

