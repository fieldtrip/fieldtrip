function matrixToNSx(varargin)

% matrixToNSx
% 
% Converts a MATLAB data matrix (channels as rows, samples as columns) into
% a NSx file.
%
% Use matrixToNSx(inputData, samplingFreq, inputUnits, savedDataPath)
% 
% All input arguments are optional.
%
%   inputData:     data matrix (channels as rows, samples as columns) into
%                  a NSx file.
%
%   samplingFreq:  Sampling frequency of the data. Currently only 500Hz,
%                  1kHz, 2kHz, 10kHz, and 30kHz sampaling frequencies are
%                  supported.
%
%   inputUnits:    The unit for the recorded data. Most data acquisition 
%                  systems save the data in units of µV. TDT units are V.
%                  The supported units are V, mV, uV or nV.
%
%   savedDataPath: The full path (excluding the extension) of the saved
%                  file.
%
%   Example 1 (Mac): 
%   matrixToNSx(myData, 30000, 'uv', '/Desktop/newConvertedFile');
%
%   In the example above, the data matrix myData recorded at 30kHz will be
%   saved on the Desktop (Mac style path) as newConvertedFile.ns5 file.
%
%   Example 2 (PC): 
%   matrixToNSx(myData, 30000, 'uv', 'C:\Data\newConvertedFile');
%
%   In the example above, the data matrix myData recorded at 30kHz will be
%   saved C:\Data folder as newConvertedFile.ns5 file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   kian@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        disp('Data is a required input.');
        return
    case 1
        disp('Sampling rate and input unit were not provided. Assuming 30kHz and µV.');
        disp('The converted file will be saved on the Desktop.');
        inputData = varargin{1};
        samplingFreq = input('What is the sampling frequency of the recorded data (e.g. ''30000'')? ');
        signalUnit = input('What is the input unit of the signal (V, mV, uV or nV)? ', 's');
        savedFilePath = input('Enter the full path for the converted file? ', 's');
    case 2
        disp('Input unit was not provided. Assuming µV.');
        disp('The converted file will be saved on the Desktop.');
        inputData = varargin{1};
        signalUnit = input('What is the input unit of the signal (V, mV, uV or nV)? ', 's');
        savedFilePath = input('Enter the full path for the converted file? ', 's');
    case 3
        disp('The converted file will be saved on the Desktop.');
        inputData = varargin{1};
        samplingFreq = varargin{2};
        signalUnit = varargin{3};
        savedFilePath = input('Enter the full path for the converted file? ', 's');
    case 4
        inputData = varargin{1};
        samplingFreq = varargin{2};
        signalUnit = varargin{3};
        savedFilePath = varargin{4};
    otherwise
        disp('Invalid number of arguments.');
        return;
end

%% Validate sampling frequency
switch samplingFreq
    case 500
        savedFileExt = '.ns1';
    case 1000
        savedFileExt = '.ns2';
    case 2000
        savedFileExt = '.ns3';
    case 10000
        savedFileExt = '.ns4';
    case 30000
        savedFileExt = '.ns5';
    otherwise
        disp('Currently only sampling rates 500Hz, 1kHz, 2kHz, 10 kHz and 30kHz are supported.');
        return;
end

%% Finding the unit of the data.
switch lower(signalUnit)
    case 'v'
        inputData = inputData * 1000000;
    case 'mv'
        inputData = inputData * 1000;
    case 'uv'
        
    case 'nv'
        inputData = inputData / 1000;
    otherwise
        disp('Invalid unit.');
        return;
end

%% Converting samplingFreq to type double
samplingFreq = double(samplingFreq);

% Warren from More lab @ Stanford
%% Determining channel count and data length
channelCount = min(size(inputData));
dataLength = max(size(inputData));

%% Creating a label of 16 character long.
newChanLabel = [num2str(samplingFreq/1000) ' kS/s'];
newChanLabel = [newChanLabel repmat(' ', [1, 16-length(newChanLabel)])];

%% Converting inputData to int16
inputData = int16(inputData);

%% Rewriting the metatags
NSx = openNSx(which('bnsx.dat'));
NSx.MetaTags.SamplingLabel = newChanLabel(1:16);
NSx.MetaTags.ChannelCount = channelCount;
NSx.MetaTags.SamplingFreq = samplingFreq;
NSx.MetaTags.TimeRes = samplingFreq;
NSx.MetaTags.ChannelID = 1:NSx.MetaTags.ChannelCount;
NSx.MetaTags.DateTime = datestr(now);
d = datevec(date);
NSx.MetaTags.DateTimeRaw = [d(1:2), weekday(date), d(3:end), 0];
NSx.MetaTags.Comment = ['This data was converted into a NSx by matrixToNSx' repmat(' ',1, 207)];
NSx.MetaTags.FileSpec = '2.3';
NSx.MetaTags.Timestamp = 0;
NSx.MetaTags.DataPoints = dataLength;
NSx.MetaTags.DataDurationSec = NSx.MetaTags.DataPoints / NSx.MetaTags.SamplingFreq;
NSx.Data = inputData;

saveNSx(NSx, [savedFilePath savedFileExt]);
