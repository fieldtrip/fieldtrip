function [hdr, participantId, data, channels, markers, impedances] = loadbvrf(hdrFile, sampleRange, flagMetadata, verbose)
% readBVRF() - Read data saved in the Brain Vision Recorder format.
%
% Usage:
%   >> [hdr, participantId, data, channels, markers, impedances] = loadbvrf(hdrFile, verbose);
%
% Inputs:
%   hdrFile      - Header file (bvrh). Data (.bvrd) and marker (bvrm) files should be in the same folder.
%   sampleRange  - [first, last] samples that define the interval to read. Enter [] to read all samples.
%   flagMetadata - Loads only the header metadata.
%   verbose      - Optional input. Set to true to display messages. 
%
% Outputs:
%   hdr             - header structure
%   participantId   - cell array of participant ids
%   data            - cell array of data matrix [channels, samples]
%   channels        - cell array of channel structure
%   markers         - cell array of marker structure
%   impedances      - cell array of impedance measurements
%
% Copyright (C) 2025 Brain Products GmbH

if nargin < 1
    error('Not enough arguments');
end

if nargin < 2
    sampleRange = [];
end

if nargin < 3
    flagMetadata = false;
end

if nargin < 4
    verbose = false;
end


if ~exist(hdrFile, "file")
    error("Input header file doesn't exist");
end

% Sanitize and read header
if verbose
    disp(['Reading: ' hdrFile]);
end
txt = fileread(hdrFile);
txt(regexp(txt, '\n')) = [];
txt(regexp(txt, ' ')) = [];
c = '﻿';
txt(txt==c) = [];
txt(txt == newline) = [];
try
    hdr = jsondecode(txt);
catch
    error("Can't parse the header")
end

if flagMetadata
    data = [];
    channels = []; 
    markers = [];
    impedances = [];
    return;
end

% Check required fields
if ~isfield(hdr, 'BVRF')
    error('Required BVRF field is missing');
end
if ~isfield(hdr, 'BVRFHeaderFile')
    error('Required BVRFHeaderFile field is missing');
end
if ~isfield(hdr, 'EEGModality')
    error('Required EEGModality field is missing');
end
if ~isfield(hdr, 'EEGModality')
    error('Required field EEGModality is missing');
end
try
    precision = hdr.EEGModality.BVRFFiles.DataFile.NumericDataType;
    precision = lower(precision);
catch
    error('Required field EEGModality.BVRFFiles.DataFile.NumericDataType is misssing')
end
try
    srate = hdr.EEGModality.DataSpecific.SamplingFrequencyInHertz;
catch
    error('Required field EEGModality.DataSpecific.SamplingFrequencyInHertz is misssing')
end
try
    channels = hdr.EEGModality.Channels;
    nChannels = length(channels);
catch 
    error('Required field EEGModality.Channels is misssing')
end


% Read data
[folder, fname, ~] = fileparts(hdrFile);
dataFile = fullfile(folder, [fname '.bvrd']);
if ~exist(dataFile, 'file')
    error('Data file is missing');
end
if verbose
    disp(['Reading: ' dataFile]);
end

if isempty(sampleRange)
    try
        fid = fopen(dataFile, 'r');
        vec = fread(fid, precision, 'ieee-le');
        fclose(fid);
        n = length(vec);
        leftOver = mod(n, nChannels);
        nSamples = (n-leftOver)/nChannels;
        data = reshape(vec(1:end-leftOver), nChannels, nSamples);
        clear vec;
    catch
        fclose(fid);
        fid = fopen(dataFile, 'r');
        n = srate*nChannels;
        data = [];
        while ~feof(fid)
            vec = fread(fid, n, precision, 'ieee-le');
            m = length(vec);
            if m < n
                leftOver = mod(m,nChannels);
                vec = reshape(vec(1:end-leftOver), nChannels, []);
            else
                vec = reshape(vec, nChannels, []);
            end
            data = [data vec];
        end
        fclose(fid);
    end

else
    if length(sampleRange) < 2
        error('Invalid sampleRange interval')
    end

    precisionMap = struct('int16', 2, 'int32', 4, 'single', 4, 'double', 8);
    if ~isfield(precisionMap, precision)
        error('EEGModality.BVRFFiles.DataFile.NumericDataType is not supported by the reader')
    end
    nBytes = precisionMap.(precision);
    fid = fopen(dataFile, 'r');
    fseek(fid, sampleRange(1) * nChannels * nBytes, 'bof');
    vec = fread(fid, sampleRange(2), precision, 'ieee-le');
    fclose(fid);
    n = length(vec);
    leftOver = mod(n, nChannels);
    nSamples = (n-leftOver)/nChannels;
    data = reshape(vec, nChannels, nSamples);
    clear vec;
end

% Post-processing
if (isstruct(channels))
    channelCell = cell(nChannels, 1);
    for k=1:nChannels
        channelCell{k} = channels(k);
    end
    channels = channelCell;
    clear channelCell;
end

for k=1:nChannels
    if isfield(channels{k}, 'Coefficients')
        warning([channels{k}.Name ' -> ' 'post-processing Coefficients is not supported']);
    end
    if isfield(channels{k}, 'ResolutionPerBit')
        data(k, :) = channels{k}.ResolutionPerBit * data(k, :);
    end
end

if isfield(channels{1}, 'ParticipantId')
    participantId = cell(nChannels, 1);
    for k=1:nChannels
        participantId{k} = channels{k}.ParticipantId;
    end
    uniqueIds = unique(participantId);
    nParticipants = length(uniqueIds);
    dataCell = cell(nParticipants, 1);
    channelsCell = cell(nParticipants, 1);
    for k=1:nParticipants
        idx = ismember(participantId, uniqueIds{k});
        dataCell{k} = data(idx, :);
        channelsCell{k} = channels(idx);
    end
    data = dataCell;
    channels = channelsCell;
    participantId = uniqueIds;
    clear dataCell channelsCell uniqueIds;
else
    participantId = {'participant-1'};
    data = {data};
    channels = {channels};
end
nParticipants = length(data);

% Read markers
markerFile = fullfile(folder, [fname '.bvrm']);
marker = struct();
markers = [];
if exist(markerFile, 'file')
    if verbose
        disp(['Reading: ' markerFile])
    end
    lines = readlines(markerFile);
    fieldNames = strsplit(lines(1), '\t', 'CollapseDelimiters',false);
    nCols = length(fieldNames);
    for fieldName = fieldNames
        marker.(fieldName) = [];
    end
    for l=2:length(lines)
        try
            content = strsplit(lines(l), '\t', 'CollapseDelimiters',false);
            if length(content) ~= nCols
                continue;
            end
            for f=1:nCols
                marker.(fieldNames(f)) = content(f);
            end
            if l>2
                markers(end+1) = marker;
            else
                markers = marker;
            end
        catch ME
            disp(ME);
            error("Invalid marker line");
        end
    end

    if nParticipants > 1
        markersCell = cell(nParticipants, 1);
        markerMask = false(nParticipants, length(markers));
        for k=1:length(markers)
            if isempty((markers(k).ParticipantId)) || contains(lower(markers(k).Comment), 'all participants')
                markerMask(:, k) = true;

            else
                idx = ismember(participantId, (markers(k).ParticipantId));
                markerMask(idx, k) = true;
            end
        end
        % Flag markers not claimed by any participant for assigment to all
        markerMask(:, all(markerMask==0)) = true;
        for k=1:nParticipants
            markersCell{k} = markers(markerMask(k, :));
        end
        markers = markersCell;
        clear markersCell;
    else
        markers = {markers};
    end
end

% Read impedance
impFile = fullfile(folder, [fname '.bvri']);
impedances = [];
if exist(impFile, 'file')
    if verbose
        disp(['Reading: ' impFile])
    end
    lines = readlines(impFile);
    if (contains(lines(1), 'ParticipantId') && nParticipants==1)
        lines(1) = [];
    elseif (contains(lines(1), 'ParticipantId') && nParticipants > 1)
        warning("Can't parse ParticipantId in the impedance file");
    end
    header = strsplit(lines(1), '\t', 'CollapseDelimiters',false);
    nCols = length(header);
    impedances = struct;
    for k=1:nCols
        impedances.(header(k)) = {};
    end
    
    for l=2:length(lines)
        try
            content = strsplit(lines(l), '\t', 'CollapseDelimiters',false);
            if length(content) ~= nCols
                continue;
            end
            for k=1:nCols
                impedances.(header(k)){end+1} = content(k);
            end
        catch ME
            disp(ME);
            error("Invalid marker line");
        end
    end 
end
end
