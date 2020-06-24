classdef KTNSPOnline

% UEAMAPFILE -- Defines a class that contains information on a UEA
% electrode array using CMP files provided by Blackrock Microsystems.
%
% To load a CMP file type 'MapName = UEAMapFile' where MapName is the name
% of the variable that will be created and will contain the map class.
%
% Example: myArrayMap = UEAMapFile;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METHODS:
%
% There are several available methods for mapfiles. Here's a list and a
% description of what they perform.
%
% Electrode2Channel(Electrode): 
%   Will return the Channel number corresponding to a passed Electrode.
%
% Channel2Electrode(Channel):   
%   Will return the Electrode number corresponding to a passed Channel.
%
% GetChannelBankID(Channel):
%   Returns the BankID for a passed Channel.
%
% GetChannelPin(Channel):
%   Returns the Pin for a passed Channel.
%
% GetChannelLabel(Channel):
%   Returns the Label for a passed Channel.
%
% GenerateChannelSubplot(Channel):
%   Returns a subplot handle for the passed Channel to be used in plotting
%   UEA-like maps.
%
% GenerateChannelSubplotNames(Channel):
%   Plots the names of channels and electrodes on the subplot corresponding
%   to the passed Channel.
%
% GetChannelColumnRow(Channel):
%   Returns the Column and Row positions of the passed number that is used
%   to plot UEA-like maps.
%
% PlotCMP:
%   Will plot a map of the CMP Map with all the electrode and channel
%   names.
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Version 1.1.0.0


%%
    properties (Hidden)
        NSPHandle
        dataCollectionFlag
        contCollection
        allSpikeData
        allLFPData
    end
    methods (Hidden)
        function obj = KTNSPOnline
            disp('Initializing NSP Connection...');
            try
                cbmex('open');
            catch
                fprintf(2,'Error connecting to NSP... Please make sure the NSP or nPlayServer is running and try again.\n');
                return;
            end
            obj.contCollection = 0;
            obj.dataCollectionFlag = 0;
            obj.dataCollectionStart;
        end
    end
    methods
        %% NSP Interface Mthods
        function closeConnection(obj)
            cbmex('close');
        end
        function openConnection(obj)
            cbmex('open');
        end
        function NSPTime = getNSPTimeSeconds(obj)
            NSPTime = cbmex('time');
        end
        function NSPTime = getNSPTimeSamples(obj)
            NSPTime = cbmex('time') * 30000;
        end
        
        %% Data Recording Methods
        function recordingStart(obj, fileName, userComment)
            if ~exist('userComment', 'var')
                userComment = 'KTNSPOnline triggered.';
            end
            if ~exist('fileName', 'var')
                fprintf(2,'Error...\n');
                disp('Filename is a required parameter.');
            else
                cbmex('fileconfig', fileName, userComment, 1);
                pause(0.1);
                cbmex('fileconfig', fileName, userComment, 1);
            end
        end
        function recordingStop(obj)
            cbmex('fileconfig', '', '', 0);
        end
        
        %% Output Methods
        function sendNSPComment(obj, comment)
            if ischar(comment)
                cbmex('comment', 0, comment);
            else
                fprintf(2,'Error...\n');
                disp('The comment should be an ASCII string.');
            end
        end
        function sendNSPMarker(obj, markerValue)
            if isnumeric(markerValue) && markerValue <= 255
                cbmex('comment', markerValue);
            else
                fprintf(2,'Error...\n');
                disp('The marker value should be a number less than 255.');
            end
        end
        function setDigOutTTLHigh(obj, digOutPortNum)
            digOutPortNum = digOutPortNum + 152;
            cbmex('digitalout', digOutPortNum, 1);
        end
        function setDigOutTTLLow(obj, digOutPortNum)
            digOutPortNum = digOutPortNum + 152;
            cbmex('digitalout', digOutPortNum, 0);
        end
        function setDigOutPulse(obj, digOutPortNum, pulseCount, frequency, pulseWidth)
            if ~exist('digOutPortNum', 'var') || ...
                    ~exist('pulseCount', 'var') || ...
                    ~exist('frequency', 'var') || ...
                    ~exist('pulseWidth', 'var')
                fprintf(2,'The function requires the following inputs:\n');
                disp('1. Digital Output (digOutPortNum): The digital out port to be triggered.');
                disp('2. Pulse Count (pulseCount):       The number of pulses to be sent out of digital out.');
                disp('3. Frequency (frequency):          The frequency of the output signal.');
                disp('4. Pulse Width (pulseWidth):       The width of the output TTL pulse.');
                disp('Example: object.setDigOutPulse(1, 10, 100, 0.01);');
            else
                digOutPortNum = digOutPortNum + 152;
                for idx = 1:pulseCount
                    cbmex('digitalout', digOutPortNum, 1);
                    pause(pulseWidth);
                    cbmex('digitalout', digOutPortNum, 0);
                    pause(1/frequency - pulseWidth);
                end
            end    
        end
        function setAnaOutPulse(obj, anaOutPortNum, pulseCount, frequency, pulseWidth)
            if ~exist('anaOutPortNum', 'var') || ...
                    ~exist('pulseCount', 'var') || ...
                    ~exist('frequency', 'var') || ...
                    ~exist('pulseWidth', 'var')
                fprintf(2,'The function requires the following inputs:\n');
                disp('1. Analog Output (anaOutPortNum): The analog out port to be triggered.');
                disp('2. Pulse Count (pulseCount):      The number of pulses to be sent out of analog out.');
                disp('3. Frequency (frequency):         The frequency of the output signal.');
                disp('4. Pulse Width (pulseWidth):      The width of the output TTL pulse in milliseconds.');
                disp('Example: object.setAnaOutPulse(1, 10, 100, 0.01);');
            else
                anaOutPortNum = anaOutPortNum + 144;
                cbmex('analogout', anaOutPortNum, 'pulses', [pulseCount 0 pulseWidth pulseWidth 4000 0 0 0], 'ms', 'mv')
            end                
        end
        
        %% Data Collection Methods
%         function status = dataCollectionIsActive(obj)
%             try
%                 activeFlag = cbmex('trialdata', 1)
%             catch
%                 status = 0;
%                 disp('Data collection is not active.');
%                 return;
%             end
%             disp('Data is currently being collected...');
%             status = 1;
%         end
        function status = dataCollectionStart(BA)
%             if ~BA.dataCollectionIsActive
                disp('Starting to collect data now...');
                cbmex('trialconfig', 1);
%             end
%             BA.dataCollectionFlag = 1;
        end
        function dataCollectionStop(obj)
%             if obj.dataCollectionIsActive
                disp('Data collection is now stopped.');
                cbmex('trialconfig', 0);
%             end
%             obj.dataCollectionFlag = 0;
        end
%         function dataCollectionSetContActive(obj)
%             disp('Continuous data collection is now active.');
%             fprintf(2, 'Do not change any channel sampling frequency while continuous data collection is active.\n');
%             obj.ContCollection = 1;
%         end
        function [spikeData, startTime, LFPData] = dataCollectionReadBuffer(obj)    
            try
                [spikeData, startTime, LFPData] = cbmex('trialdata', 1);
            catch
                fprintf(2,'Error...\n');
                disp('Data collection is not active.');
                disp('Use dataCollectionStart method to start locating.');
            end
%             if obj.contCollection
%                 numChannelsRead = size(LFPData,1);
%                 if isempty(obj.allLFPData)
%                     obj.allLFPData = LFPData;
%                 else
%                     for idx = 1:numChannelsRead
%                         obj.allLFPData{idx, 3} = [obj.allLFPData{idx, 3}; LFPData{idx, 3}];
%                     end
%                 end
%             end
        end
%         function data = getDataLFPContinuous(obj)
%             data = obj.allLFPData;
%         end
        
        %% Detection Methods
        function spikedChannels = detectChannelsFired(obj, units)
            if ~exist('units', 'var')
                units = 1:5;
            end
            spikeData = dataCollectionReadBuffer(obj);
            spikeData(:,1) = [];
            spikedChannels = ~cellfun(@isempty, spikeData(:,units));
        end
        function fireFlag = detectChannelUnitFiredAny(obj, channels, units)
            if ~exist('units', 'var')
                units = 1:5;
            end
            spikedChannels = obj.detectChannelsFired(units);
            fireFlag = any(any(spikedChannels(channels, units)));
        end
        function fireFlag = detectChannelUnitFiredAll(obj, channels, units)
            if ~exist('units', 'var')
                units = 1:5;
            end
            spikedChannels = obj.detectChannelsFired(units);
            fireFlag = all(all(spikedChannels(channels, units)));
        end
        function fireFlag = detectWhenChannelUnitFiredAny(obj, channels, units, timeWindow)
            fireFlag = 0;
            if ~exist('units', 'var')
                units = 1:5;
            end
            if ~exist('timeWindow', 'var')
                timeWindow = 0.02;
            end
            while ~fireFlag
                pause(timeWindow);
                fireFlag = detectChannelUnitFiredAny(obj, channels, units);
            end
        end
%         function fireFlag = detectWhenChannelUnitFiredAll(obj, channels, units, timeWindow)
%             fireFlag = 0;
%             if ~exist('units', 'var')
%                 units = 1:5;
%             end
%             if ~exist('timeWindow', 'var')
%                 timeWindow = 0.02;
%             end
%             while ~fireFlag
%                 pause(timeWindow);
%                 fireFlag = detectChannelUnitFiredAll(obj, channels, units);
%             end
%         end
        function [bitValues, bitTimestamps] = detectDigInWord(obj)
            readData = obj.dataCollectionReadBuffer;
            bitValues = readData{151,3};
            bitTimestamps = readData{151,2};
        end
        function [bitValues, bitTimestamps] = detectDigInBinary(obj)
            readData = obj.dataCollectionReadBuffer;
            bitValues = dec2bin(readData{151,3});
            bitTimestamps = readData{151,2};
        end
% % % %         function fireFlag = detectDigInBit(obj, bitDetect)
% % % %             fireFlag = 0;
% % % %             if ~exist('bitDetect', 'var')
% % % %                 disp('Please specify a bit that needs to be detected.');
% % % %                 disp('Example: detectDigInBit(3)');
% % % %                 return;
% % % %             end
% % % %             spikedChannels = obj.detectChannelsFired([1:5]);
% % % %             if (bitDetect < 0) || (bitDetect > 15)
% % % %                 disp('The bit value can be between 0 and 15, inclusive.');
% % % %             else
% % % %                 fireFlag = any(any(spikedChannels(136, 1:5)));
% % % %             end
% % % %         end
        function fireFlag = detectWhenDigInWord(obj, wordValue, timeWindow)
            fireFlag = 0;
            if ~exist('wordValue', 'var')
                fprintf(2,'Error...\n');
                disp('Please specify a word value that needs to be detected.');
                disp('Example: detectDigInBit(128)');
                return;
            end
            if ~exist('timeWindow', 'var')
                timeWindow = 0.02;
            end
            while ~fireFlag
                pause(timeWindow);
                readWordValue = detectDigInWord(obj);
                if wordValue == readWordValue
                    fireFlag = 1;
                end
            end
        end
        function fireFlag = detectWhenDigInBit(obj, bitsValue, timeWindow)
            fireFlag = 0;
            if ~exist('bitValue', 'var')
                fprintf(2,'Error...\n');
                disp('Please specify a bit pattern that needs to be detected.');
                disp('Example: detectDigInBit(13)');
                return;
            end
            if ~exist('timeWindow', 'var')
                timeWindow = 0.02;
            end
            while ~fireFlag
                pause(timeWindow);
                bitsValue = detectDigInBit(obj, bitDetect);
                if strcmpi(bitsValue, readWordValue)
                    fireFlag = 1;
                end
            end
        end
        
        %% Get Configuration Values
        function configValue = getAllChannelLabels(obj)
            configStream = cbmex('chanlabel');
            configValue = configStream(:,1);
        end
        function configValue = getChannelLabel(obj, channel)
            configStream = cbmex('chanlabel', channel);
            configValue = configStream{1};
        end
        function isFlag = isChannelSpikeExtractionEnabled(obj, channel)
            configStream = cbmex('chanlabel', channel);
            isFlag = configStream{2};
        end
        function isFlag = isChannelUnitEnabled(obj, channel, unit)
            configStream = cbmex('chanlabel', channel);
            isFlag = configStream{unit};
        end
        function ampRejRange = getAmplitudeRejectionRange(obj, channel)
            configStream = cbmex('config', channel);
            ampRejRange = [configStream{12,2}, configStream{13,2}];               
        end
        function setAmplitudeRejectionRange(obj, channel, newValue)
            if ~exist('newValue', 'var')
                fprintf(2,'Error...\n');
                disp('newValue is a required argument.');
                return;
            end
            if length(newValue)<2 || length(newValue)>2
                fprintf(2,'Error...\n');
                disp('The passed argument to this function should have two values: lower and higher amplitude rejection values, respectively.');
                return;
            end
            if newValue(1)>=newValue(2)
                fprintf(2,'Error...\n');
                disp('The low amplitude rejection value should be smaller than the high value.');
                return;
            end
            configStream = cbmex('config', channel, 'amplrejneg', newValue(1));
            configStream = cbmex('config', channel, 'amplrejpos', newValue(2));          
        end       
        function spikeThresh = getSpikeThresholdValue(obj, channel)
            configStream = cbmex('config', channel);
            spikeThresh = ceil(configStream{10, 2}/4);
        end
        function setSpikeThresholdValue(obj, channel, newValue)
            if ~exist('newValue', 'var')
                fprintf(2,'Error...\n');
                disp('newValue is a required argument.');
                return;
            else
                configStream = cbmex('config', channel, 'spkthrlevel', newValue * 4);
            end
        end
        function refChan = getDigitalReferenceChannel(obj, channel)
            configStream = cbmex('config', channel);
            refChan = configStream{14,2};
        end
        function refChan = setDigitalReferenceChannel(obj, channel, newValue)
            if ~exist('newValue', 'var')
                fprintf(2,'Error...\n');
                disp('newValue is a required argument.');
                return;
            else
                configStream = cbmex('config', channel, 'refelecchan', newValue);
            end
        end
        
        %% Data Analysis Methods
        function onlineSpectogram(obj, frequencyRange)
            if ~exist('frequencyRange', 'var')
                fprintf(2,'Error...\n');
                disp('You need to define a frequency range for this function.');
                disp('Example: onlineSpectrogram([0.3, 7500]) to show the spectrogram in the range of 0.3 Hz to 7500 Hz.');
                return;
            end
            frequencyRange = linspace(frequencyRange(1), frequencyRange(end), 150);
            % All durations are in second
            sampleCollectionDuration = 0.02;
            refreshRate = 0.05;
            % Setting up a default frequency range, if none is specified
            if ~exist('frequencyRange', 'var')
                disp('Frequency range was not specified. A default range of 0.3:1:7500 will be used.');
                frequencyRange = 0.3:1:7500;
            end
            % Setting up the figure
            procFigure = figure;
            set(procFigure, 'Name', 'Close this figure to stop'); 
            % Setting up collection variable
            timeDisplay = tic;
            timeCollection = tic;
            flagCollection = 1;
            % Collection
            dataCollectionStart(obj);
            while ishandle(procFigure)
                pause(0.001);
                if flagCollection
                    timeCollectionET = toc(timeCollection);
                    if timeCollectionET > sampleCollectionDuration
                        [spikeData, startTime, contData] = dataCollectionReadBuffer(obj);
                        nGraphs = size(contData,1);
                        if ishandle(procFigure)
                            for idx = 1:nGraphs
                                fs0 = contData{idx, 2};
                                data = contData{idx, 3};
                                collectSize = min(size(data), sampleCollectionDuration * fs0);
                                x = data(1:collectSize);
                                if isempty(frequencyRange)
                                    [psd, f] = periodogram(double(x),[],'onesided',512,fs0); 
                                else
                                    [psd, f] = periodogram(double(x),[],frequencyRange,fs0);
                                end
                                subplot(nGraphs,1,idx,'Parent',procFigure);
                                psdLog = 10*log10(psd);
                                plot(f, psdLog, 'b');
                                title(sprintf('fs = %d t = %f',fs0, startTime));
                                xlabel('Frequency (Hz)'); xlim([frequencyRange(1), frequencyRange(end)]);
                                ylabel('Magnitude (dB)');
                            end
                            drawnow;
                        end
                        flagCollection = 0;
                    end
                end
                timeCollectionET = toc(timeDisplay);
                if timeCollectionET >= refreshRate;
                    timeDisplay = tic;
                    timeCollection = tic;
                    flagCollection = 1;
                end
            end
            dataCollectionStop(obj);
        end       
    end
end