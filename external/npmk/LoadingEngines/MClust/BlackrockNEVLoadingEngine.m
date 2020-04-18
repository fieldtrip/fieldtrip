function [timeStamps, waveForms] = BlackrockNEVLoadingEngine(varargin)

%% 
% Blackrock Microsystems NEV loading engine for MClust 3.0.
%
% Kian Torab
% support@blackrockmicro.com
% Blackrock Microsystems
% 
% Version 1.3.0.0
%
% Version 1.4 by Saman Hagh Gooie - 13 October 2016
%       - Modifed to work with the latest version of the MClust v4.4


if nargin == 1
    fn = varargin{1};
elseif nargin == 2
    if strcmp(varargin{1}, 'get')
        switch (varargin{2})
            case 'ChannelValidity'
                timeStamps = [true true true true]; return;
            case 'ExpectedExtension'
                timeStamps = '.nev'; return;
            case 'UseFileDialog'
                timeStamps = true; return;
            case 'filenames'
                timeStamps = false; return;
            otherwise
                error('Unknown get condition.');
        end
    else
        error('2 argins requires "get" as the first argument.');
    end
elseif nargin == 3
    fn = varargin{1};
    recordToGet = varargin{2};
    recordUnits = varargin{3};
end


if ~exist(fn, 'file')
    disp('File does not exist.');
    timeStamps = false;
    waveForms = false;
    return;
end
clc
%% Reading the entire data if recordUnits is not defined
if ~exist('recordUnits', 'var')
    NEV = loadNEV(fn, 'read');
    %Determining what channels are in the data file
    readElectrodes = unique(NEV.Data.Spikes.Electrode);    
    
    %Saving all timestamps from the file into the timeStamp output variable
    timeStamps = NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==readElectrodes(1));
    timeStamps = timeStampsToTimestampSeconds(timeStamps,NEV.MetaTags.SampleRes);
    %Saving all saveforms from the file into the waveForm output variable
    waveForms = initializeWaveform(length(timeStamps));
   
    for idx = 1:length(readElectrodes)
        waveForms(:,idx,:) = NEV.Data.Spikes.Waveform(1:32,NEV.Data.Spikes.Electrode == readElectrodes(idx))';
    end
else
%% Reading the file depending on what recordUnits is defining
    switch recordUnits
        case 1 %Returning timestamp and waveforms for given timestamps
            %Loading the NEV file
            NEV = loadNEV(fn, 'read');
            %Saving all timestamps in the timeStamps output variable
            timeStamps = NEV.Data.Spikes.TimeStamp;
            %Converting seconds to NEV timestamps
            recordToGet = timeStampsSecondsToTimestamps(recordToGet,NEV.MetaTags.SampleRes);
            %Calculating the index of timestamps to get
            for idx = 1:length(recordToGet)
                IDXtoGet(idx) = find(timeStamps == recordToGet(idx),1);
            end
            %Saving all timestamps requested by the timestamps
            timeStamps = timeStamps(IDXtoGet);
            timeStamps = timeStampsToTimestampSeconds(timeStamps,NEV.MetaTags.SampleRes);
            %Calculating the number of timestamps
            numOfTimestamps  = length(timeStamps);
            %Preallocating a 3D waveforms matrix
            waveForms = initializeWaveform(numOfTimestamps);
            %Saving all waveforms in the waveForms output variable
            waveForms(:,1,1:32) = NEV.Data.Spikes.Waveform(1:32,IDXtoGet)';
        case 2 %Returning timestamp and waveforms for given records
            %Loading the NEV file
            NEV = loadNEV(fn, 'read');
            %Saving the timestamps for the requested records
            timeStamps = NEV.Data.Spikes.TimeStamp(1,recordToGet);
            timeStamps = timeStampsToTimestampSeconds(timeStamps,NEV.MetaTags.SampleRes);
            %Calculating the number of timestamps
            numOfTimestamps = length(timeStamps);
            %Saving the waveforms for the requested records
            waveForms = initializeWaveform(numOfTimestamps);
            %Determining what channels are in the data file
            readElectrodes = unique(NEV.Data.Spikes.Electrode);    
            for idx = 1:length(readElectrodes)
                temp = NEV.Data.Spikes.Waveform(1:32,NEV.Data.Spikes.Electrode == readElectrodes(idx))';
                waveForms(:,idx,:) = temp(recordToGet,:);
            end
        case 3 %Returning timestamp and waveforms for range of given timestamps
            %Loading the NEV file
            NEV = loadNEV(fn, 'read');
            %Converting timestamps in seconds to timestamps in samples
            recordToGet = timeStampsSecondsToTimestamps(recordToGet,NEV.MetaTags.SampleRes);
            %Parsing out the timestamp range
            tInitial = recordToGet(1);
            tEnd = recordToGet(2);
            %Saving all timestamps in the timeStamps output variable
            timeStamps = NEV.Data.Spikes.TimeStamp;
            %Finding timestamps that fall within the range
            IDX = find(timeStamps <= tEnd & timeStamps >= tInitial);
            %Saving all timestamps that fall within the range
            timeStamps = NEV.Data.Spikes.TimeStamp(1,IDX);
            timeStamps = timeStampsToTimestampSeconds(timeStamps,NEV.MetaTags.SampleRes);
            %Calculating the number of timestamps
            numOfTimestamps  = length(timeStamps);
            %Preallocating a 3D waveforms matrix
            waveForms = initializeWaveform(numOfTimestamps);
            %Determining what channels are in the data file
            readElectrodes = unique(NEV.Data.Spikes.Electrode);    
            for idx = 1:length(readElectrodes)
                temp = NEV.Data.Spikes.Waveform(1:32,NEV.Data.Spikes.Electrode == readElectrodes(idx))';
                waveForms(:,idx,:) = temp(IDX,:);
            end
        case 4 %Returning timestamp and waveforms for the range of given records
            %Loading the NEV file
            NEV = loadNEV(fn, 'read');
            %Creating a range out of the record list
            recordToGet = [recordToGet(1):recordToGet(end)];
            %Saving the timestamps for the requested records
            timeStamps = NEV.Data.Spikes.TimeStamp(1,recordToGet);
            timeStamps = timeStampsToTimestampSeconds(timeStamps,NEV.MetaTags.SampleRes);
            %Calculating the number of timestamps
            numOfTimestamps = length(timeStamps);
            %Saving the waveforms for the requested records
            waveForms = initializeWaveform(numOfTimestamps);
            % Loading the waveforms into a MClust waveform variable
            readElectrodes = unique(NEV.Data.Spikes.Electrode);
            for idx = 1:length(readElectrodes)
                waveForms(:,idx,:) = NEV.Data.Spikes.Waveform(1:32,NEV.Data.Spikes.Electrode == readElectrodes(idx))';
            end
            %Saving the waveforms for the requested records
            waveForms = waveForms(recordToGet, :, :);
            
        case 5 %Returning the number of records in the file
            %Loading the NEV file
            NEV = loadNEV(fn);
            % Returning the number of spikes
            timeStamps = size(NEV.Data.Spikes.TimeStamp, 2)/length(unique(NEV.Data.Spikes.Electrode));
            waveForms = [];
        otherwise
            disp('Invalid input for recordUnit');
            disp('Values from 1 to 5 are valid only.');
    end
end

timeStamps = timeStamps';

%% Uses openNEV (a part of NPMK package: http://bit.ly/brnpmkkit) to open the NEV file
function NEV = loadNEV(fn, input)

if ~exist('input', 'var')
    NEV = openNEV(fn, 'nomat', 'nosave');
    
elseif strcmpi(input, 'noread')
    NEV = openNEV(fn, 'nomat', 'nosave');
    
elseif ~strcmpi(input, 'read')
    disp('Invalid input.');
    
else
    NEV = openNEV(fn, 'nomat', 'nosave');
end

%% Initializes the waveform variable depending on the number of spikes in the file
function waveForms = initializeWaveform(numOfTimestamps)

    waveForms = zeros(numOfTimestamps, 4, 32);
    
%% Converts the timestamp samples to timestamp values in seconds
function timeStampsSec = timeStampsToTimestampSeconds(timeStamps, samplingRate)

    timeStampsSec = double(timeStamps) / double(samplingRate);

%% Converts the timestamp seconds to timestamp values in samples
function timeStamps = timeStampsSecondsToTimestamps(timeStampsSeconds, samplingRate)

    timeStamps = int32(timeStampsSeconds * double(samplingRate));
        
        
