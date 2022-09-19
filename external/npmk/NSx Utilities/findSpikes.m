function Spikes = findSpikes(NSx, varargin)
% findSpikes
%
% Searches NSx data structure for spikes by thresholding. The output is
% compatible with NEV.Data.Spikes data structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use OUTPUT = findSpikes(NSx, 'threshold', 'preThreshold', 'spikeLength', 'channels', 'duration', 'filter')
%   
%    NSx:          The data structure holding the NSx structure
%
%    NOTE: All following input arguments are optional. Input arguments may 
%          be in any order.
%
%   'threshold':    The threshold value to apply to the file.
%                   DEFAULT: -65 uV
%
%   'preThreshold': The number of pre-threshold crossing samples to save.
%                   DEFAULT: 10
%
%   'spikeLength':  Length of each extracted spike in number of samples.
%                   DEFAULT: 48
%
%   'channels':     The channels to perform spike extraction on.
%                   DEFAULT: will perform spike extraction on all channels
%
%   'duration':     The duration of file to perform spike extraction on.
%                   DEFAULT: will perform extraction on the entire file
%
%   'filter':       The filter cut-off frequency applied before spike
%                   extraction.
%                   DEFAULT: no filter will be applied
%
%   OUTPUT:      Contains the NEV.Spikes structure.
%
%   USAGE EXAMPLE: 
%   
%   Spikes = findSpikes(NSx, 'channels', 1:3, 'duration', 1:10000, 'threshold', -65);
%   
%   In the above example the NSx analog data is searched for spikes on
%   channels 1 through 3 for samples between 1 and 10000. Threshold is
%   set to -65 uV for all channels.
%
%   Spikes = findSpikes(NSx, 'channels', 1:3, 'threshold', [-65 -100 -85]);
%   In the above example the NSx analog data is searched for spikes among
%   channels 1 to 3. Thresholds for channel 1 is set to -65 uV, for channel
%   2 is set to -100 uV, and for channel 3 is set to -85 uV.
%
%   Original Author: Ehsan Azar
%
%   Contributors: 
%   Kian Torab, Blackrock Microsystems, ktorab@blackrockmicro.com
%
%   Version 1.0.2.0
%

Spikes = struct('TimeStamp', [],'Electrode', [], 'Unit', [],'Waveform', [], 'findSpikesVer', []);

%% Validating the input arguments. Exit with error message if error occurs.
spikelen = 48;
threshold = -65;
preThreshold = 10;
filterCorner = 250;
next = '';
for i=1:length(varargin)
    inputArgument = varargin{i};
    if (strcmpi(next, 'threshold'))
        next = '';
        threshold = inputArgument;
    elseif (strcmpi(next, 'preThreshold'))
        next = '';
        preThreshold = inputArgument;
    elseif (strcmpi(next, 'channels'))
        next = '';
        channels = inputArgument;
    elseif (strcmpi(next, 'duration'))
        next = '';
        duration = inputArgument;
    elseif (strcmpi(next, 'spikeLength'))
        next = '';
        spikelen = inputArgument;
    elseif (strcmpi(next, 'filter'))
        next = '';
        filterCorner = inputArgument;
    elseif strcmpi(inputArgument, 'threshold')
        next = 'threshold';
    elseif strcmpi(inputArgument, 'preThreshold')
        next = 'preThreshold';        
    elseif strcmpi(inputArgument, 'channels')
        next = 'channels';        
    elseif strcmpi(inputArgument, 'duration')
        next = 'duration';        
    elseif strcmpi(inputArgument, 'spikeLength')
        next = 'spikeLength';        
    elseif strcmpi(inputArgument, 'filter')
        next = 'filter'; 
    end
end
clear next;

%% Give all input arguments a default value. All input argumens are
%  optional.
if ~exist('channels', 'var'); channels = (1:size(NSx.Data, 1))'; end
if ~exist('duration', 'var'); duration = (1:size(NSx.Data, 2))'; end
if ~exist('filter', 'var');   filter = [];  end

%% Apply filter
Data = NSx.Data(channels, duration).';

if ~isempty(filter)
    [b, a] = butter(4, filterCorner/15000, 'high');
    Data = filter(b, a, Data);
end

%% Threshold
if isfield(NSx, 'ElectrodesInfo');
    if (isscalar(threshold) && all([NSx.ElectrodesInfo(channels).MaxAnalogValue] == NSx.ElectrodesInfo(1).MaxAnalogValue) && ...
            all([NSx.ElectrodesInfo(channels).MaxDigiValue] == NSx.ElectrodesInfo(1).MaxDigiValue))
        threshold = (threshold * double(NSx.ElectrodesInfo(1).MaxDigiValue) / ...
            double(NSx.ElectrodesInfo(1).MaxAnalogValue)); % scalar threshold
    else
        threshold = (threshold .* double([NSx.ElectrodesInfo(channels).MaxDigiValue]) ./ ...
            double([NSx.ElectrodesInfo(channels).MaxAnalogValue])); % vector threshold
    end
    threshold = cast(threshold, class(Data));
end

[row, col] = find((diff(bsxfun(@minus, Data, threshold) > 0) < 0).');
clear threshold;
Spikes.TimeStamp = col' - preThreshold;
clear col;
Spikes.Electrode = row';
clear row;

Spikes.Electrode(Spikes.TimeStamp < 1) = [];
Spikes.TimeStamp(Spikes.TimeStamp < 1) = [];
Spikes.Electrode(Spikes.TimeStamp + spikelen > size(NSx.Data, 2)) = [];
Spikes.TimeStamp(Spikes.TimeStamp + spikelen > size(NSx.Data, 2)) = [];

%% Lockout violation removal
% this makes sure all spikes are at least one spike length apart
ii = 1;
while (ii < length(Spikes.TimeStamp))
    idx = find((Spikes.Electrode((ii+1):end) == Spikes.Electrode(ii)) & ...
        (Spikes.TimeStamp((ii+1):end) - Spikes.TimeStamp(ii) >= spikelen), 1);
    if (isempty(idx))
        idx = find(Spikes.Electrode((ii+1):end) == Spikes.Electrode(ii));
    else
        idx = find(Spikes.Electrode((ii+1):(ii+idx-1)) == Spikes.Electrode(ii));
    end
    if (~isempty(idx))
        Spikes.Electrode(ii + idx) = [];
        Spikes.TimeStamp(ii + idx) = [];
    end
    ii = ii + 1;
end
clear idx;

%% Spike extraction
Spikes.Waveform = zeros(spikelen, length(Spikes.TimeStamp));
for ii = 1:length(Spikes.TimeStamp)
    ts = Spikes.TimeStamp(ii):(Spikes.TimeStamp(ii) + spikelen - 1);
    Spikes.Waveform(:, ii) = Data(ts, Spikes.Electrode(ii));
end
% convert index to actual channel number
Spikes.Electrode = channels(Spikes.Electrode);


