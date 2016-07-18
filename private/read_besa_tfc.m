function [ChannelLabels, Time, Frequency, Data, Info] = read_besa_tfc(FILENAME)

% READ_BESA_TFC imports data from a BESA *.tfc file
%
% Use as
%   [DataType, ConditionName, Channels, Time, Frequency, Data] = read_besa_tfc(FILENAME)
%
% This reads data from the BESA Time-Frequency-Coherence output data file
% FILENAME and returns the following data:
%   ConditionName: name of analyzed condition
%   ChannelLabels: character array of channel labels
%   Time: array of sampled time instants
%   Frequency: array of sampled frequencies
%   Data: 3D data matrix with indices (channel,time,frequency)
%   Info: Struct containing additional information:
%       DataType: type of the exported data
%       ConditionName: name of analyzed condition
%       NumbeOfTrials: Number of trials on which the data is based
%       StatisticsCorrection: Type of statistics correction for multiple testing
%       EvokedSignalSubtraction: Type of evoked signal subtraction 

% Copyright (C) 2005, Vladimir Litvak
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

fp = fopen(FILENAME);
   
VersionNumber = fscanf(fp,'VersionNumber=%s ',1);
DataType = fscanf(fp,'DataType=%s ',1);
ConditionName = fscanf(fp,'ConditionName=%s ',1);
try
    if ConditionName == 'Condition'
        ConditionName = [ConditionName,' ',fscanf(fp,'%s ',1)];
    end
catch
end
NumberTrials = fscanf(fp,'NumberTrials=%i ',1);
NumberTimeSamples = fscanf(fp,'NumberTimeSamples=%i ');
TimeStartInMS = fscanf(fp,'TimeStartInMS=%f ',1);
TimeIntervalInMS = fscanf(fp,'IntervalInMS=%f ',1);
NumberFrequencies = fscanf(fp,'NumberFrequencies=%i ');
FreqStartInHZ = fscanf(fp,'FreqStartInHz=%f ',1);
FreqIntervalInHZ = fscanf(fp,'FreqIntervalInHz=%f ',1);
NumberChannels = fscanf(fp,'NumberChannels=%i ');

% New file versions (BESA 5.0.8 and higher) include more information in the .tfc file header; skip that
vers=0;                 % new tfc file format
try
    StatisticsCorrection = fscanf(fp,'StatisticsCorrection=%s ',1);
    EvokedSignalSubtraction = fscanf(fp,'EvokedSignalSubtraction=%s',1);
catch
    vers=1;             % old tfc file format
end

% Handle possible future extensions of the tfc file header
i=1;
while i<1000
    a = fscanf(fp,'%c',1);
    if strcmp(a,sprintf('\n'))
        i=1000;
    end
    i=i+1;
end
    
% Generate return values
% FIXME the following statement does not work for MATLAB 7.2 on XP (see mail from Stephan Bickel)
Time = [TimeStartInMS:TimeIntervalInMS:(NumberTimeSamples-1)*TimeIntervalInMS+TimeStartInMS];
Frequency = [FreqStartInHZ:FreqIntervalInHZ:(NumberFrequencies-1)*FreqIntervalInHZ+FreqStartInHZ];
if vers == 1
    Info = struct('DataType',{DataType},'ConditionName',{ConditionName},'NumberOfTrials',{NumberTrials});
else
    Info = struct('DataType',{DataType},'ConditionName',{ConditionName},'NumberOfTrials',{NumberTrials},...
        'StatisticsCorrection',{StatisticsCorrection},'EvokedSignalSubtraction',{EvokedSignalSubtraction});
end

if isempty(findstr(DataType,'COH'))
    for Channel=1:NumberChannels
        ChannelLabels(Channel) = cellstr(fscanf(fp,'%s ',1));
    end
else
    for Channel=1:NumberChannels
        ChannelLabels(Channel) = cellstr(fscanf(fp,'%s ',3));
    end    
end

ChannelLabels=char(ChannelLabels);

Data = zeros(NumberChannels,NumberTimeSamples,NumberFrequencies);
for Channel=1:NumberChannels
    Data(Channel,:,:) = fscanf(fp,'%f',[NumberTimeSamples,NumberFrequencies]);
end
fclose(fp);
