% pop_loadbvrf() - Import BrainVision Recording Format (BVRF) dataset 
%                  into EEGLAB. Supports multi-participant recordings, 
%                  partial channel import, sample-range import, and optional 
%                  import of markers, impedances, and metadata.
%
% Usage:
%   >> [ALLEEG, com] = pop_loadbvrf(hdrPath, hdrFileName);
%   >> [ALLEEG, com] = pop_loadbvrf(hdrPath, hdrFileName, 'key', value, ...);
%
% Graphic interface:
%   >> [ALLEEG, com] = pop_loadbvrf;  
%      (opens GUI to select header and import options)
%
% Required Inputs:
%   hdrPath     - Path to the *.bvrh header file.
%   hdrFileName - Name of the header file (string, '*.bvrh').
%
% Optional Key/Value Inputs:
%   'sampleInterval'       - [firstSample lastSample] to read. Empty [] reads full dataset.
%
%   'channelIndx'          - Vector of channel indices to import. If empty [], all channels are imported.
%                             Channel order corresponds to the "Channels" section of the BVRF header
%
%   'flagImportMarkers'    - true/false. Import markers from *.bvrm file.(Default: True)
%
%   'flagImportImpedances' - true/false. Import impedance information from *.bvri file? 
%                            This file is CONDITIONAL and only present if
%                            impedances were recorded. (Default: False)
%
%   'participantId'          - If dataset includes multiple participants, 
%                              you can select one by provindg its ID. Empty [] reads all.
%                             (Deault: [])
%
%   'usePoly'                - true/false Use polynomial
%                              information in sensors (if exist) to convert sensor 
%                              output to physical quantity. (Default: true)
%
% Outputs:
%   ALLEEG - Cell array of EEGLAB EEG structures (one per selected participant)
%   com    - Command string to reproduce the import
%
% Example:
%   >> ALLEEG = pop_loadbvrf('/data/subj1/', 'recording.bvrh', ...
%                            'sampleInterval', [0 100000], ...
%                            'channelIndx', 1:32, ...
%                            'flagImportMarkers', true);
%
% Author: Ramon Martinez-Cancino, Brain Products GmbH, 2025
%         Alejandro Ojeda, Brain Products GmbH, 2025
%
% Copyright (C) 2025 Brain Products GmbH
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [ALLEEG, com] = pop_loadbvrf(hdrPath, hdrFileName,varargin) 

com = '';
EEG = []; ALLEEG = {};

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            opt.(options{i}) = options{i+1};
        end
    else
        opt = [];
    end
catch
    disp('pop_pac() error: calling convention {''key'', value, ... } error'); return;
end
try opt.participantId;           catch, opt.participantId         = [];     end
try opt.sampleInterval;          catch, opt.sampleInterval        =  [];    end
try opt.channelIndx;             catch, opt.channelIndx           =  [];    end
try opt.flagImportMarkers;       catch, opt.flagImportMarkers     =  true;  end
try opt.flagImportImpedances;    catch, opt.flagImportImpedances  =  false; end
try opt.usePoly;                 catch, opt.usePoly               = true;   end


% Pick Header file if not provided
if nargin < 2
    [hdrFileName, hdrPath] = uigetfile2({'*.bvrh'}, 'Select bvrh-file - pop_loadbvrf()');
    if hdrFileName(1) == 0, return; end
end

% Load Header file
disp('pop_loadbvrf(): reading header file');
hdrFileFullPath = fullfile(hdrPath, hdrFileName);
hdr = loadbvrf(hdrFileFullPath, [], 1);

% Launch UI if needed
if nargin < 2
    %-- info structure for UI
    info = struct();
    info.bvrfVersion   = hdr.BVRF.Version;                                           % BVRF Version
    info.dataType      = hdr.EEGModality.BVRFFiles.DataFile.NumericDataType;         % Data type                                                 
    info.srate         = hdr.EEGModality.DataSpecific.SamplingFrequencyInHertz;      % Hz sampling rate

    % Determine samples
    [filepath, filename]= fileparts(hdrFileFullPath);
    dataFile = fullfile(filepath, [filename, '.bvrd']);
    if ~exist(dataFile, 'file')
        error('Data file is missing');
    end

    % 1. Number of channels
    numChannels = numel(hdr.EEGModality.Channels);

    % 2. Bytes per sample based on NumericDataType
    bytesPerSample = struct('Int16',2,'Int32',4,'Single',4,'Double',8);
    bps = bytesPerSample.(info.dataType);

    % 3. DataFile size
    fileInfo = dir(dataFile);
    fileBytes = fileInfo.bytes;

    % 4. Compute number of samples
    info.samples = floor(fileBytes / (numChannels * bps));

    % Determine participants
    if isfield(hdr, 'Participants')
        % Multiple participants explicitly listed
        Participants = normalizeObjectArray(hdr.Participants,'Participants,');
        for isubj =1: numel(Participants)
            participantIds{isubj} = Participants{isubj}.Id;
        end
        % participantIds = {hdr.Participants.Id};
        nParticipants = numel(participantIds);
    else
        % Single participant (Participants section omitted)
        participantIds = {'participant-1'};
        nParticipants = 1;
    end
    info.participantNames = participantIds;

    % Extract channels section
    channels = normalizeObjectArray(hdr.EEGModality.Channels, 'EEGModality.Channels');

    % Assign channels to participants
    if nParticipants == 1
        % All channels belong to participant 1
        chanCounts = numel(channels);
    else
        % Count per participant using ParticipantId in channels
        chanCounts = zeros(1, nParticipants);
        uniqueIds = cellfun(@(p) p.Id, hdr.Participants, 'UniformOutput', false);

        for ichan =1:numel(channels)
            ChannelsParticipantID{ichan} = safe_get(channels{ichan},'ParticipantId', '');
        end
        isUnassigned = cellfun('isempty', ChannelsParticipantID);

        for p = 1:nParticipants
            idx = ismember(ChannelsParticipantID, uniqueIds{p}) | isUnassigned;
            chanCounts(p) = nnz(idx);
        end
    end
    info.participantChans = chanCounts;

    % Additional data present
    % Detect markers

    info.hasMarkers = isfield(hdr.EEGModality.BVRFFiles, 'MarkerFile');

    % Detect impedance file
    info.hasImpedances = isfield(hdr.EEGModality.BVRFFiles, 'ImpedanceFile');
    % --
    cfg = bvrf_reader_gui(info); % Launch UI
    if isempty(cfg), return, end

    opt.sampleInterval = cfg.sampleInterval;
    opt.channelIndx = cfg.channelIndx;
    opt.flagImportMarkers = cfg.flagImportMarkers;
    opt.flagImportImpedances = cfg.flagImportImpedances';
    opt.participantId = cfg.participantId;
    opt.usePoly = cfg.usePoly;
end

%% HDR to EEG structure matching (From here on , the code is independent of the calling way)
 [~, ALLEEGtmp] = eeg_loadbvrf(hdrPath, hdrFileName,'sampleInterval', opt.sampleInterval,...
                                                         'channelIndx', opt.channelIndx, ...
                                                         'flagImportMarkers', opt.flagImportMarkers,...
                                                         'flagImportImpedances', opt.flagImportImpedances,...
                                                         'participantId', opt.participantId,...
                                                         'usePoly', opt.usePoly,...
                                                         'verbose',true);

for isubj = 1:numel(ALLEEGtmp)
    try
        ALLEEG{isubj} = eeg_checkset(ALLEEGtmp{isubj}); %#ok<AGROW>
    catch
    end
end

if nargout == 2
    com = sprintf('ALLEEG = pop_loadbvrf(''%s'', ''%s'',''participantID'', %s, ''sampleInterval'', %s, ''channelIndx'',%s, ''flagImportMarkers'',  %s, ''flagImportImpedances'', %s, ''usePoly'', %s);', ...
           hdrPath, hdrFileName, opt.participantId, mat2str(opt.sampleInterval), mat2str(opt.channelIndx), mat2str(opt.flagImportMarkers),mat2str(opt.flagImportImpedances), mat2str(opt.usePoly));
end
