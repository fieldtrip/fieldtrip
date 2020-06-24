function NEV = rethresholdNSx(varargin)

% rethresholdNSx
% 
% Opens a NEV file and rethresholds it. It outputs a NEV structure into
% MATLAB. This version currently does not save a new NEV file.
% 
% All input arguments are optional. Input arguments must be in the
% following order:
%
%   NSx:          The data structure holding the channel information
%
%   Threshold:    The threshold value in µV. This can be a positive or
%                 negative value.
%                 DEFAULT: Will automatically set threshold to -65 µV.
%
%   FilterFreq:   The high-pass filter high-pass corner frequency. For
%                 example, if the filter is set to 250 then the signal will
%                 first get high-pass filtered at 250 Hz before getting
%                 rethrholded.
%                 DEFAULT: Will automatically high-pass filter signal at
%                 250 Hz.
%
%   Example: 
%   
%   NEV = rethresholdNSx(NSx, -85, 150);
%
%   In the example above, the NSx structure will be high-pass filtered at 
%   150 Hz and then rethresholded at -85 µV. All new waveformes will be
%   saved in a NEV structure and will get output to MATLAB.
%
%   NEV = rethresholdNSx;
%
%   In the example above, the user will be propted to choose a NSx file
%   then the signal will get high-pass filtered at 150 Hz and then 
%   rethresholded at -85 µV. All new waveformes will be saved in a NEV 
%   structure and will get output to MATLAB.
%
%   Kian Torab
%   Blackrock Microsystems
%   kian@blackrockmicro.com
%
%   Version 1.0.0.0
%

if nargin > 3
    disp('Invalid number of arguments. See ''help rethresholdNSx'' for more information.');
elseif nargin == 3
    NSx = varargin{1};
    try
        NEV = openNEV([NSx.MetaTags.FilePath '/' NSx.MetaTags.Filename(1:end-3) 'nev'], 'read');
    catch
        disp('Cannot open the corresponding NEV file.');
        disp('The NEV file is required. Please place it in the same folder.');
        return;
    end
    threshold = varargin{2};
    hpFreq = varargin{3};
elseif nargin == 2
    NSx = varargin{1};
    try
        NEV = openNEV([NSx.MetaTags.FilePath '/' NSx.MetaTags.Filename(1:end-3) 'nev'], 'read');
    catch
        disp('Cannot open the corresponding NEV file.');
        disp('The NEV file is required. Please place it in the same folder.');
        return;
    end
    threshold = varargin{2};
    hpFreq = 250;
elseif nargin == 1
    NSx = varargin{1};
    try
        NEV = openNEV([NSx.MetaTags.FilePath '/' NSx.MetaTags.Filename(1:end-3) 'nev'], 'read');
    catch
        disp('Cannot open the corresponding NEV file.');
        disp('The NEV file is required. Please place it in the same folder.');
        return;
    end
    threshold = -85;
    hpFreq = 250;
elseif nargin == 0
    [fileName pathName] = getFile;
    NSx = openNSx([pathName fileName], 'read');
    try
        NEV = openNEV([pathName fileName(1:end-3) 'nev'], 'read');
    catch
        disp('Cannot open the corresponding NEV file.');
        disp('The NEV file is required. Please place it in the same folder.');
        return;
    end
    threshold = -50;
    hpFreq = 250;
end

%% Setting stuff up
newSpikes.TimeStamps = 0;
newSpikes.Electrode = 0;
newSpikes.Unit = 0;
newSpikes.Waveform = zeros(48,1);
firstTimeFlag = 1;

digitizationFactor = 1000/NEV.ElectrodesInfo(1).DigitalFactor;
wavelengthLength = length(NEV.Data.Spikes.Waveform(:,1));
availableChannels = double(NSx.MetaTags.ChannelID)';

%% Highpass filter continuous data



%% Finding the spikes

for channelIDX = 1:length(availableChannels)
    if threshold < 0
        thresholdTimestamps = find(diff(NSx.Data(channelIDX,:) < threshold*digitizationFactor) == 1);
    else
        thresholdTimestamps = find(diff(NSx.Data(channelIDX,:) > threshold*digitizationFactor) == 1);
    end
    idx = 1;
    while idx < length(thresholdTimestamps)
        refractoryIDX = find(abs(thresholdTimestamps(idx) - thresholdTimestamps) < wavelengthLength);
        thresholdTimestamps(refractoryIDX(2:end)) = [];
        idx = idx + 1;
    end
    while length(NSx.Data(channelIDX, :)) - thresholdTimestamps(end) < 38
        thresholdTimestamps(end) = [];
    end
    newSpikes.TimeStamps(end+1:end+length(thresholdTimestamps)) = thresholdTimestamps;
    newSpikes.Electrode(end+1:end+length(thresholdTimestamps)) = repmat(NSx.MetaTags.ChannelID(channelIDX), 1, length(thresholdTimestamps));
    newSpikes.Unit(end+1:end+length(thresholdTimestamps)) = zeros(1, length(thresholdTimestamps));
    if firstTimeFlag
        newSpikes.TimeStamps(1) = [];
        newSpikes.Electrode(1) = [];
        newSpikes.Unit(1) = [];
    end
    waveformIndices = bsxfun(@plus, thresholdTimestamps', -10:37)';
    waveformFlat = NSx.Data(channelIDX, waveformIndices(:));
    newSpikes.Waveform(:,end+1:end+length(thresholdTimestamps)) = reshape(waveformFlat, 48, length(waveformFlat)/48);
    if firstTimeFlag
        newSpikes.Waveform(:,1) = [];
        firstTimeFlag = 0;
    end
    clear waveformFlat;
end

[~, timestampIndex] = sort(newSpikes.TimeStamps);

NEV.Data.Spikes.TimeStamp = newSpikes.TimeStamps(timestampIndex);
NEV.Data.Spikes.Electrode = newSpikes.Electrode(timestampIndex);
NEV.Data.Spikes.Unit = newSpikes.Unit(timestampIndex);
NEV.Data.Spikes.Waveform = newSpikes.Waveform(:,timestampIndex);