function tfc = readBESAtfc(filename)

% readBESAtfc reads all information from a *.tfc file
%
% Use as
%   tfc = readBESATFC(filename)
%
% The output is a structure containing the following fields:
%   ChannelLabels: character array of channel labels
%   Time: array of sampled time instants
%   Frequency: array of sampled frequencies
%   Data: 3D data matrix with indices (channel,time,frequency)
%   DataType: type of the exported data
%   ConditionName: name of analyzed condition
%   NumberOfTrials: Number of trials on which the data is based
%   StatisticsCorrection: Type of statistics correction for multiple testing
%   EvokedSignalSubtraction: Type of evoked signal subtraction
%
% Last modified April 25, 2006 Karsten Hoechstetter
% Last modified April 26, 2006 Robert Oostenveld

if isempty(findstr(filename,'.'))
  filename = [filename,'.tfc'];
end
fp = fopen(filename);

VersionNumber = fscanf(fp,'VersionNumber=%s ');
tfc.DataType = fscanf(fp,'DataType=%s ');
tfc.ConditionName = fscanf(fp,'ConditionName=%s ');

% If the user has not specified a condition name, BESA defaults to
% "Condition 1" etc. This is the only instance where the condition
% name can have a blank
if strcmp(tfc.ConditionName,'Condition')
  tfc.ConditionName = [tfc.ConditionName,' ',num2str(fscanf(fp,'%i '))];
end
temp = fscanf(fp,'NumberTrials=%f NumberTimeSamples=%f TimeStartInMS=%f IntervalInMS=%f NumberFrequencies=%f FreqStartInHz=%f FreqIntervalInHz=%f NumberChannels=%i ' );
tfc.NumberOfTrials = temp(1);
NumberTimeSamples = temp(2);
TimeStartInMS = temp(3);
TimeIntervalInMS = temp(4);
NumberFrequencies = temp(5);
FreqStartInHZ = temp(6);
FreqIntervalInHZ = temp(7);
NumberChannels = temp(8);

% New file versions (BESA 5.0.8 and higher) include more information in the .tfc file header; skip that
vers=0;                 % new tfc file format
try
  tfc.StatisticsCorrection = fscanf(fp,'StatisticsCorrection=%s ');
  tfc.EvokedSignalSubtraction = fscanf(fp,'EvokedSignalSubtraction=%s');
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
tfc.Time = TimeStartInMS:TimeIntervalInMS:(NumberTimeSamples-1)*TimeIntervalInMS+TimeStartInMS;
tfc.Frequency = FreqStartInHZ:FreqIntervalInHZ:(NumberFrequencies-1)*FreqIntervalInHZ+FreqStartInHZ;

if isempty(findstr(tfc.DataType,'COH'))
  for Channel=1:NumberChannels
    tfc.ChannelLabels(Channel) = cellstr(fscanf(fp,'%s ',1));
  end
else
  for Channel=1:NumberChannels
    tfc.ChannelLabels(Channel) = cellstr(fscanf(fp,'%s ',3));
  end
end

tfc.ChannelLabels=char(tfc.ChannelLabels);

tfc.Data = zeros(NumberChannels,NumberTimeSamples,NumberFrequencies);
for Channel=1:NumberChannels
  tfc.Data(Channel,:,:) = fscanf(fp,'%f',[NumberTimeSamples,NumberFrequencies]);
end
fclose(fp);
