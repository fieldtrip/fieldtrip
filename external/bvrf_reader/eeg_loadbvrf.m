function [hdr, ALLEEG] = eeg_loadbvrf(hdrPath, hdrFileName, varargin)
% eeg_loadbvrf() - Read data saved in the Brain Vision Recorder format (BVRF)
%                  and create one EEGLAB EEG struct per participant.
%
% Usage:
%   >> [hdr, ALLEEG] = eeg_loadbvrf(hdrPath, hdrFileName, 'key', value, ...);
%
% Inputs:
%   hdrPath      - Header path
%   hdrFileName  - Header file name (.bvrh). Data (.bvrd), marker (.bvrm),
%                  and impedance (.bvri) files must be in the same folder.
%
% Optional 'key', value pairs:
%   'sampleInterval'       - [first, last] sample indices to read (0-based).
%                            [] (default) = read all samples.
%   'channelIndx'          - Indices of channels to import (applied per participant).
%   'flagImportMarkers'    - true/false (default: true). Import markers → EEG.event.
%   'flagImportImpedances' - true/false (default: true). Read impedance file
%                            and store in EEG.etc.impedances.
%   'flagMetadata'         - true/false (default: false). If true: only read
%                            header; ALLEEG is returned empty.
%   'verbose'              - true/false (default: false). Print progress.
%   'participantId'        - string. If provided and there are multiple
%                            participants, only the participant whose
%                            Channels(k).ParticipantId matches this value
%                            is loaded into ALLEEG.
%   'usePoly'              - true/false (default: true). Use polynomial
%                            information in sensors (if exist) to convert sensor
%                            output to physical quantity
%
% Outputs:
%   hdr        - Parsed BVRF header (struct from jsondecode).
%   ALLEEG     - 1 x N cell array of EEGLAB EEG structs (one per participant,
%                or just one if 'participantId' is specified).
%
% Notes:
%   - Participant splitting is done via Channels(k).ParticipantId only.
%   - Markers are mapped to EEGLAB events as:
%         EEG.event(k).type    = MarkerType
%         EEG.event(k).latency = Sample
%         EEG.event(k).comment = Comment
%
% Author: Alejandro Ojeda, Brain Products GmbH, 2025
%         Ramon Martinez-Cancino, Brain Products GmbH, 2025
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

ALLEEG = {};

%% ---------------------------------------------------------------------
%   Parse inputs
% ----------------------------------------------------------------------
if nargin < 2
    error('eeg_loadbvrf:NotEnoughArguments', 'Not enough arguments.');
end

% Options
opt = struct;
try
    options = varargin;
    if ~isempty(options)
        for i = 1:2:numel(options)
            opt.(options{i}) = options{i+1};
        end
    end
catch
    disp('eeg_loadbvrf() error: calling convention {''key'', value, ... } error');
    return;
end

% Defaults
try opt.sampleInterval;          catch, opt.sampleInterval        = [];     end
try opt.channelIndx;             catch, opt.channelIndx           = [];     end
try opt.flagImportMarkers;       catch, opt.flagImportMarkers     = true;   end
try opt.flagImportImpedances;    catch, opt.flagImportImpedances  = false;  end
try opt.flagMetadata;            catch, opt.flagMetadata          = false;  end
try opt.verbose;                 catch, opt.verbose               = false;  end
try opt.participantId;           catch, opt.participantId         = [];     end
try opt.usePoly;                 catch, opt.usePoly               = true;   end

%% ---------------------------------------------------------------------
%   File checks
% ----------------------------------------------------------------------
[~, fname, ~] = fileparts(hdrFileName);
hdrFile    = fullfile(hdrPath, hdrFileName);
dataFile   = fullfile(hdrPath, [fname '.bvrd']);
markerFile = fullfile(hdrPath, [fname '.bvrm']);
impFile    = fullfile(hdrPath, [fname '.bvri']);

if ~exist(hdrFile,  "file"), error("Input header file doesn't exist"); end
if ~exist(dataFile, 'file') && ~opt.flagMetadata
    error('Data file is missing');
end

if opt.flagImportMarkers && ~exist(markerFile, 'file')
    warning('Marker file is missing; markers will not be imported.');
end
if opt.flagImportImpedances && ~exist(impFile, 'file')
    warning('Impedance file is missing; impedances will not be imported.');
end

%% ---------------------------------------------------------------------
%   Read and parse header
% ----------------------------------------------------------------------
if opt.verbose
    disp(['Reading header: ' hdrFile]);
end
txt = fileread(hdrFile);
txt(regexp(txt, '\n')) = [];
txt(regexp(txt, ' '))  = [];
c = '﻿';
txt(txt == c)       = [];
txt(txt == newline) = [];

try
    hdr = jsondecode(txt);
catch
    error("Can't parse the header");
end

% Required fields
if ~isfield(hdr, 'BVRF')
    error('Required BVRF field is missing');
end
if ~isfield(hdr, 'BVRFHeaderFile')
    error('Required BVRFHeaderFile field is missing');
end
if ~isfield(hdr, 'EEGModality')
    error('Required EEGModality field is missing');
end

try
    precision = hdr.EEGModality.BVRFFiles.DataFile.NumericDataType;
    precision = lower(precision);
catch
    error('Required field EEGModality.BVRFFiles.DataFile.NumericDataType is missing');
end

try
    srate = hdr.EEGModality.DataSpecific.SamplingFrequencyInHertz;
catch
    error('Required field EEGModality.DataSpecific.SamplingFrequencyInHertz is missing');
end

try
    channelsHdr = hdr.EEGModality.Channels;
    nChannels   = numel(channelsHdr);
catch
    error('Required field EEGModality.Channels is missing');
end

% If only metadata requested, stop here.
if opt.flagMetadata
    if opt.verbose
        disp('flagMetadata is true: returning header only');
    end
    ALLEEG = {};
    return;
end

%% ---------------------------------------------------------------------
%   Read data (.bvrd)
% ----------------------------------------------------------------------
if opt.verbose
    disp(['Reading data: ' dataFile]);
end

if isempty(opt.sampleInterval)
    % Read all samples (try fast path)
    try
        fid = fopen(dataFile, 'r');
        vec = fread(fid, precision, 'ieee-le');
        fclose(fid);
        n        = length(vec);
        leftOver = mod(n, nChannels);
        nSamples = (n - leftOver) / nChannels;
        data     = reshape(vec(1:end-leftOver), nChannels, nSamples);
        clear vec;
    catch
        % Fallback: read in blocks
        fclose(fid);
        fid = fopen(dataFile, 'r');
        nPerBlock = srate * nChannels;
        data = [];
        while ~feof(fid)
            vec = fread(fid, nPerBlock, precision, 'ieee-le');
            m   = length(vec);
            if m == 0
                break;
            end
            if m < nPerBlock
                leftOver = mod(m, nChannels);
                vec = reshape(vec(1:end-leftOver), nChannels, []);
            else
                vec = reshape(vec, nChannels, []);
            end
            data = [data vec]; %#ok<AGROW>
        end
        fclose(fid);
    end
else
    % Read only requested sample interval
    if numel(opt.sampleInterval) < 2
        error('Invalid sampleInterval; must be [first last].');
    end

    precisionMap = struct('int16', 2, 'int32', 4, 'single', 4, 'double', 8);
    if ~isfield(precisionMap, precision)
        error('EEGModality.BVRFFiles.DataFile.NumericDataType is not supported by the reader');
    end
    nBytes = precisionMap.(precision);

    fid = fopen(dataFile, 'r');
    fseek(fid, opt.sampleInterval(1) * nChannels * nBytes, 'bof');
    vec = fread(fid, opt.sampleInterval(2) * nChannels, precision, 'ieee-le');
    fclose(fid);

    n        = length(vec);
    leftOver = mod(n, nChannels);
    nSamples = (n - leftOver) / nChannels;
    data     = reshape(vec(1:end-leftOver), nChannels, nSamples);
    clear vec;
end

%% ---------------------------------------------------------------------
%   Post-processing: normalize channels + scaling
% ----------------------------------------------------------------------
channels = normalizeObjectArray(hdr.EEGModality.Channels, 'Channels');

for k = 1:nChannels
    ResolutionPerBit = safe_get(channels{k}, 'ResolutionPerBit', 1);
    % Extract coefficients (if any)
    [hasPoly, num, denom] = bvrfGetChannelCoeffs(channels{k});

    if hasPoly && opt.usePoly
        % MQ(NV) = polynomial(NV, ResolutionPerBit, Coeficients)
        data(k, :) = evalRationalPolyAscending(double(data(k, :)), num, denom, ResolutionPerBit);
        if opt.verbose
            fprintf('Coefficients applied to Channel: %d \n', k);
        end
    else
        % MQ = NV * RE
        data(k, :) = double(data(k, :)) .* ResolutionPerBit;
    end
end

%% ---------------------------------------------------------------------
%   Split data and channels by ParticipantId
% ----------------------------------------------------------------------
if isfield(hdr, 'Participants')
    uniqueIds = cellfun(@(p) p.Id, hdr.Participants, 'UniformOutput', false);
    nParticipants = numel(uniqueIds);

    for ichan =1:numel(channels)
        participantId{ichan} = safe_get(channels{ichan},'ParticipantId', '');
    end

    dataCell     = cell(nParticipants, 1);
    channelsCell = cell(nParticipants, 1);
    isUnassigned = cellfun('isempty', participantId);

    if nnz(isUnassigned)~=0 && opt.verbose
        fprintf('Unassigned channel(s) have been detected: %d.  Proceeding to append to each participant ... \n', nnz(isUnassigned));
    end

    for p = 1:nParticipants
        idx              = ismember(participantId, uniqueIds{p}) | isUnassigned;
        dataCell{p}      = data(idx, :);
        channelsCell{p}  = channels(idx);
    end

    data         = dataCell;
    channels     = channelsCell;
    participantId = uniqueIds;

    clear dataCell channelsCell uniqueIds;
else
    % Single participant case
    participantId = {'participant-1'};
    data          = {data};
    channels      = {channels};
    nParticipants = 1;
end

%% ---------------------------------------------------------------------
%   Optional: select one participant by participantId
% ----------------------------------------------------------------------
if ~isempty(opt.participantId) && nParticipants>1
    targetId = char(string(opt.participantId));
    idxKeep  = strcmp(participantId, targetId);

    if ~any(idxKeep)
        error('eeg_loadbvrf:ParticipantNotFound', ...
              'No participant with ParticipantId = "%s" found in header.', targetId);
    end

    participantId = participantId(idxKeep);
    data          = data(idxKeep);
    channels      = channels(idxKeep);
    nParticipants = numel(data); % likely 1, but allow >1 if header uses same ID multiple times

    if opt.verbose
        fprintf('Selected participantId = %s (N=%d).\n', targetId, nParticipants);
    end
end

%% Apply channel index selection per participant (if requested)
if ~isempty(opt.channelIndx)
    for p = 1:nParticipants
        chanIdx = opt.channelIndx;
        if max(chanIdx) > size(data{p},1)
            error('Requested channelIndx exceeds number of channels for participant %d.', p);
        end
        data{p}     = data{p}(chanIdx, :);
        channels{p} = channels{p}(chanIdx);
    end
end

%% ---------------------------------------------------------------------
%   Read markers (.bvrm) and split per participant
% ----------------------------------------------------------------------
markers = {};
if opt.flagImportMarkers && exist(markerFile, 'file')
    if opt.verbose
        disp(['Reading markers: ' markerFile]);
    end
    lines = readlines(markerFile);
    if isempty(lines)
        markers = cell(nParticipants,1);
    else
        fieldNames = strsplit(lines(1), '\t', 'CollapseDelimiters', false);
        nCols      = numel(fieldNames);

        markerTemplate = struct();
        for f = 1:nCols
            markerTemplate.(fieldNames(f)) = [];
        end

        markersAll = [];
        for l = 2:numel(lines)
            try
                content = strsplit(lines(l), '\t', 'CollapseDelimiters', false);
                if numel(content) ~= nCols
                    continue;
                end
                marker = markerTemplate;
                for f = 1:nCols
                    marker.(fieldNames(f)) = content(f);
                end
                if isempty(markersAll)
                    markersAll = marker;
                else
                    markersAll(end+1) = marker; %#ok<AGROW>
                end
            catch ME
                disp(ME);
                error('Invalid marker line in file.');
            end
        end

        if isempty(markersAll)
            markers = cell(nParticipants,1);
        else
            % Assign markers to participants via ParticipantId
            markersCell = cell(nParticipants, 1);
            markerMask  = false(nParticipants, numel(markersAll));

            for k = 1:numel(markersAll)
                mk = markersAll(k);

                if ~isfield(mk, 'ParticipantId') || isempty(mk.ParticipantId) || ...
                        contains(lower(string(mk.Comment)), 'all participants')
                    markerMask(:, k) = true;
                else
                    idx = ismember(participantId, mk.ParticipantId);
                    markerMask(idx, k) = true;
                end
            end

            % Any marker not claimed -> assign to all
            markerMask(:, all(markerMask == 0, 1)) = true;

            for p = 1:nParticipants
                markersCell{p} = markersAll(markerMask(p, :));
            end
            markers = markersCell;
        end
    end
elseif opt.flagImportMarkers
    % File missing but flagImportMarkers true -> keep empty cell
    markers = cell(nParticipants, 1);
end

%% ---------------------------------------------------------------------
%   Read impedances (.bvri)
% ----------------------------------------------------------------------
impedances = [];
if opt.flagImportImpedances && exist(impFile, 'file')
    if opt.verbose
        disp(['Reading impedances: ' impFile]);
    end
    lines = readlines(impFile);
    if isempty(lines)
        impedances = [];
    else
        % Handle ParticipantId in first line if needed
        if contains(lines(1), 'ParticipantId') && nParticipants == 1
            lines(1) = [];
        elseif contains(lines(1), 'ParticipantId') && nParticipants > 1
            warning("Can't parse ParticipantId in the impedance file; attaching raw table.");
        end

        header = strsplit(lines(1), '\t', 'CollapseDelimiters', false);
        nCols  = numel(header);
        impedances = struct;
        for k = 1:nCols
            impedances.(header(k)) = {};
        end

        for l = 2:numel(lines)
            try
                content = strsplit(lines(l), '\t', 'CollapseDelimiters', false);
                if numel(content) ~= nCols
                    continue;
                end
                for k = 1:nCols
                    impedances.(header(k)){end+1} = content(k);
                end
            catch ME
                disp(ME);
                error('Invalid impedance line in file.');
            end
        end
    end
end

%% ---------------------------------------------------------------------
%   Build EEG structs per participant
% ----------------------------------------------------------------------
ALLEEG = cell(1, nParticipants);

for p = 1:nParticipants
    EEG = emptyset;

    % --- core data/time fields ---
    EEG.data   = double(data{p});
    [EEG.nbchan, EEG.pnts] = size(EEG.data);
    EEG.srate  = double(srate);
    EEG.trials = 1;
    EEG.xmin   = 0;
    EEG.xmax   = EEG.pnts / EEG.srate;
    EEG.comments = [ 'Original file: ' hdrFile ];

    % --- subject / condition / setname ---
    subjId = char(string(participantId{p}));
    EEG.subject = subjId;

    taskName = '';
    if isfield(hdr, 'Task') && isfield(hdr.Task, 'Name') && ~isempty(hdr.Task.Name)
        taskName = char(string(hdr.Task.Name));
    end

    if ~isempty(taskName)
        EEG.condition = taskName;
        EEG.setname   = sprintf('%s_%s', taskName, subjId);
    else
        EEG.condition = '';
        EEG.setname   = sprintf('BVRF_%s', subjId);
    end

    % --- chanlocs and reference ---
    EEG.chanlocs = build_chanlocs_from_channels(channels{p}, hdr.EEGModality);
    EEG.ref = compute_ref(channels{p});

    % Store full BVRF reference information for later reconstruction
    EEG.etc.bvrf = struct();

    % 1) Store the named reference definitions (if provided in the header)
    if isfield(hdr.EEGModality, 'References')
        EEG.etc.bvrf.references = hdr.EEGModality.References;
    else
        EEG.etc.bvrf.references = [];
    end

    % --- markers -> EEG.event ---
    EEG.event = [];
    if opt.flagImportMarkers && ~isempty(markers) && numel(markers) >= p && ~isempty(markers{p})
        EEG.event = markers_to_eeg_events(markers{p}, EEG.srate);
    end

    % --- impedances ---
    if opt.flagImportImpedances && ~isempty(impedances)
        EEG.etc.impedances = impedances;
    end

    % Minimal bookkeeping
    EEG.etc.participant_id = subjId;
    EEG.etc.bvrf_header    = hdr;  % remove if you want lighter structs

    ALLEEG{p} = EEG;
end
end

%% =====================================================================
%   Helper: build_chanlocs_from_channels
% =====================================================================
function chanlocs = build_chanlocs_from_channels(channelsCell,eegMod)

chStruct = channelsCell;
nCh      = numel(chStruct);

chanlocs = struct('labels', cell(1,nCh), ...
    'type',   cell(1,nCh), ...
    'X',      cell(1,nCh), ...
    'Y',      cell(1,nCh), ...
    'Z',      cell(1,nCh));

% Normalize electrodes and sensors from header
electrodes = normalizeObjectArray(eegMod.Electrodes, 'Electrodes');
sensors    = normalizeObjectArray(eegMod.Electrodes, 'Sensors');

for k = 1:nCh
    ch = chStruct{k};
    chanlocs(k).labels = safe_get(ch, 'Name', sprintf('Ch%d', k));
    chanlocs(k).type   = safe_get(ch, 'Type', '');

    [X, Y, Z] = deal([]);
    if isfield(ch, 'Composition') && isfield(ch.Composition, 'Plus')
        plus = ch.Composition.Plus;
        name = safe_get(plus, 'Name', '');
        sec  = safe_get(plus, 'Section', 'Electrodes');
        [X, Y, Z] = lookup_xyz(name, sec, electrodes, sensors);
    end
    chanlocs(k).X = X;
    chanlocs(k).Y = Y;
    chanlocs(k).Z = Z;
end
end

%% =====================================================================
%   Helper: compute_ref
% =====================================================================
function refStr = compute_ref(channelsCell)
% Reference logic (EEGLAB-compliant):
%
% - Channels.Composition.Minus is OPTIONAL in BVRF.
% - “common” in EEGLAB means **all channels share the same reference**.
%
% Rules:
%   1) If NO channel has a Minus field → EEG.ref = 'unknown'
%   2) If SOME channels have Minus and others don't → EEG.ref = 'unknown'
%   3) If all channels have Minus but they differ → EEG.ref = 'unknown'
%   4) If all channels have identical Minus:
%        - If Minus.Name is non-empty → EEG.ref = that name (e.g. 'Cz', 'REF')
%        - If Minus.Name is empty     → EEG.ref = 'common' (shared but unnamed)
%
% Note:
% - EEG.ref is metadata only; EEGLAB does not use it during rereferencing.
% - Full BVRF Minus definitions is preserved in EEG.etc for accuracy.

chStruct = channelsCell;
nCh      = numel(chStruct);
if nCh == 0
    refStr = 'unknown';
    return;
end

minusInfo = cell(nCh,1);
hasMinus  = false(1,nCh);

for k = 1:nCh
    ch = chStruct{k};
    if isfield(ch,'Composition') && isfield(ch.Composition,'Minus')
        m = ch.Composition.Minus;
        minusInfo{k} = struct( ...
            'Name',    safe_get(m,'Name',''), ...
            'Section', safe_get(m,'Section','') );
        hasMinus(k) = true;
    else
        minusInfo{k} = struct('Name','', 'Section','');
    end
end

% --- If some channels have Minus and others don't → mixed, unknown
if ~all(hasMinus) && any(hasMinus)
    refStr = 'unknown';
    return;
end

% --- If no channels have Minus → unknown
if ~any(hasMinus)
    refStr = 'unknown';
    return;
end

% --- All have Minus: check if identical
same = true;
for k = 2:nCh
    if ~isequal(minusInfo{1}, minusInfo{k})
        same = false;
        break;
    end
end

if ~same
    refStr = 'unknown';
    return;
end

% --- All identical Minus
if ~isempty(minusInfo{1}.Name)
    refStr = char(string(minusInfo{1}.Name));
else
    refStr = 'common';
end

end

%% =====================================================================
%   Helper: markers_to_eeg_events
% =====================================================================
function events = markers_to_eeg_events(markerStructArray, srate)

if isempty(markerStructArray)
    events = struct('type', {}, 'latency', {}, 'comment', {});
    return;
end

nEv    = numel(markerStructArray);
events = struct('type', cell(1,nEv), ...
    'latency', cell(1,nEv), ...
    'comment', cell(1,nEv),...
    'bv_type', cell(1,nEv),...
    'bv_StartEndId',cell(1,nEv), ...
    'bv_channel', cell(1, nEv));

for k = 1:nEv
    mk = markerStructArray(k);

    % BV type
    if isfield(mk, 'Type') && ~isempty(mk.Type)
        bvtypeStr = char(string(mk.Type));
    else
        fprintf(2,'Field ''Type''  in the Marker structure is <strong>Required </strong>, but is not provided \nWe are Replacing Marker Type with generic name: ''Marker''\n')
        bvtypeStr = 'Marker';
    end

    % General type
    bvcode = safe_get(mk, 'Code', []);
    bvvalue = safe_get(mk, 'Value', []);

    if bvcode == '', bvcode = []; end
    if bvvalue == '', bvvalue = []; end

    if ~isempty(bvcode) && ~isempty(bvvalue)
        typeStr = strcat( bvtypeStr, char('/'), char(bvcode), char(num2str(bvvalue)));

    elseif ~isempty(bvcode) && isempty(bvvalue)
        typeStr = strcat( bvtypeStr, char('/'), char(bvcode));

    elseif isempty(bvcode) && ~isempty(bvvalue)
        typeStr = strcat( bvtypeStr, char('/'), char(num2str(bvvalue)));

    else
        typeStr = bvtypeStr;
    end

    % latency: prefer Sample, fallback to SampleIndex, fallback to Time [s]
    latency = NaN;
    if isfield(mk, 'Sample') && ~isempty(mk.Sample)
        latency = str2double(string(mk.Sample));

    elseif isfield(mk, 'SampleIndex') && ~isempty(mk.SampleIndex)
        latency = str2double(string(mk.SampleIndex));

    elseif isfield(mk, 'Time') && ~isempty(mk.Time)
        t = str2double(string(mk.Time));
        latency = 1 + round(t * srate); % seconds -> sample index
    end

    if isnan(latency) || ~isfinite(latency)
        continue; % skip invalid
    end

    %% Output Fields
    % Standard EEGLAB  fields
    events(k).type           = typeStr;
    events(k).latency        = latency + 1; % Adding one to account for the zero indexing of the original value.
    events(k).comment        = safe_get(mk, 'Comment', '');

    % BV Specific fields
    events(k).bv_type        = bvtypeStr;
    events(k).bv_code        = safe_get(mk, 'Code', '');
    events(k).bv_value       = safe_get(mk, 'Value', '');
    events(k).bv_StartEndId  = safe_get(mk, 'StartEndId', '');
    events(k).bv_channel     = safe_get(mk, 'Channel', '');

end

% Remove any incomplete events (where latency is NaN)
validMask = arrayfun(@(e) ~isempty(e.latency) && isfinite(e.latency), events);
events    = events(validMask);

end


%% =====================================================================
%   Helper: Get Channels coeficients
% =====================================================================
function [hasPoly, num, denom] = bvrfGetChannelCoeffs(chan)

hasPoly = false;
num   = [];
denom = [];

if ~isfield(chan, 'Coefficients') || isempty(chan.Coefficients)
    return;
end

if ~isfield(chan.Coefficients, 'Num') || isempty(chan.Coefficients.Num)
    return;
end

num = chan.Coefficients.Num(:)';

if isfield(chan.Coefficients, 'Denom') && ~isempty(chan.Coefficients.Denom)
    denom = chan.Coefficients.Denom(:)';
else
    denom = 1;   % default denominator
end

% Default identity polynomial (spec default): Num = [0 1], Denom = [1]
defaultNum = [0 1];
defaultDen = 1;

isIdentity = isequal(num, defaultNum) && isequal(denom, defaultDen) ;
hasPoly = ~isIdentity;
end

%% =====================================================================
%   Generic helpers
% =====================================================================
function [X, Y, Z] = lookup_xyz(name, section, electrodes, sensors)
X = []; Y = []; Z = [];
if isempty(name)
    return;
end

switch lower(section)
    case 'electrodes'
        arr = electrodes;
    case 'sensors'
        arr = sensors;
    otherwise
        arr = electrodes;
end

if isempty(arr)
    return;
end

names = arrayfun(@(x) arr{x}.Name, 1:numel(arr), 'UniformOutput', false);
idx   = find(strcmp(names, name), 1);
if ~isempty(idx)
    coords = safe_get(arr{idx}, 'Coordinates', []);
    if numel(coords) >= 3
        X = coords(1); Y = coords(2); Z = coords(3);
    end
end
end

function EEG = emptyset()
EEG.setname     = '';
EEG.filename    = '';
EEG.filepath    = '';
EEG.subject     = '';
EEG.group       = '';
EEG.condition   = '';
EEG.session     = [];
EEG.comments    = '';
EEG.nbchan      = 0;
EEG.trials      = 0;
EEG.pnts        = 0;
EEG.srate       = 1;
EEG.xmin        = 0;
EEG.xmax        = 0;
EEG.times       = [];
EEG.data        = [];
EEG.icaact      = [];
EEG.icawinv     = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.chanlocs    = [];
EEG.urchanlocs  = [];
EEG.chaninfo    = [];
EEG.ref         = [];
EEG.event       = [];
EEG.urevent     = [];
EEG.eventdescription = {};
EEG.epoch       = [];
EEG.epochdescription = {};
EEG.reject      = [];
EEG.stats       = [];
EEG.specdata    = [];
EEG.specicaact  = [];
EEG.splinefile  = '';
EEG.icasplinefile = '';
EEG.dipfit      = [];
EEG.history     = '';
EEG.saved       = 'no';
EEG.etc         = [];
end
