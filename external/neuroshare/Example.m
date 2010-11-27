function Example()

% Prompt for the correct DLL
disp(' ');  % Blank line
DLLName = input('DLL Name: ', 's');

% Load the appropriate DLL
[nsresult] = ns_SetLibrary(DLLName);
if (nsresult ~= 0)
    disp('DLL was not found!');
    return
end

% Find out the data file from user
disp(' ');  % Blank line
filename = input('Data file: ', 's');


% Load data file and display some info about the file
% Open data file
[nsresult, hfile] = ns_OpenFile(filename);
if (nsresult ~= 0)
    disp('Data file did not open!');
    return
end
clear filename;

% Get file information
[nsresult, FileInfo] = ns_GetFileInfo(hfile);
% Gives you EntityCount, TimeStampResolution and TimeSpan
if (nsresult ~= 0)
    disp('Data file information did not load!');
    return
end
   
% Define some variables needed for firing rates
stepsize = 0.02; % seconds
stepsize1 = 0.1; % seconds
if FileInfo.TimeSpan > 150  % Limit the timespan shown in the graphs to 150 seconds
    totaltime = 150;
else
    totaltime = FileInfo.TimeSpan; % seconds
end
time = 0 : stepsize : totaltime;   % Initialize time axis for gaussian plot

% Build catalogue of entities
[nsresult, EntityInfo] = ns_GetEntityInfo(hfile, [1 : 1 : FileInfo.EntityCount]);

NeuralList = find([EntityInfo.EntityType] == 4);    % List of EntityIDs needed to retrieve the information and data
SegmentList = find([EntityInfo.EntityType] == 3);
AnalogList = find([EntityInfo.EntityType] == 2);
EventList = find([EntityInfo.EntityType] == 1);

% How many of a particular entity do we have
cNeural = length(NeuralList);       
cSegment = length(SegmentList);
cAnalog = length(AnalogList);
cEvent = length(EventList);

clear FileInfo;

if (cNeural == 0)
    disp('No neural events available!');
end

if (cSegment == 0)
    disp('No segment entities available!');
    return;     % It does not make sense to continue in this particular analysis
                % if there are no segment entities.
end

if (cAnalog == 0)
    disp('No analog entities available!');
end

if (cEvent == 0)
    disp('No event entities available!');
end

% Have user pick a channels or channels for further analysis
% Show the user how many channels are available
disp(' ');
disp(['There are ' num2str(cSegment) ' channels.']);
disp(' ');

channel = input('Which data channels would you like to display? (e.g. 1 or [1 2 3])');
if (channel > cSegment)  % Have to check that the selected channel actually exists
    disp('Channel does not exist');
    return
end
clear cNeural cSegment;

% Throw away entity infos we don't need to save memory
EntityInfo = rmfield(EntityInfo, 'EntityType');

%
% Load the waveform data and do the analysis
%

% These three functions are for demonstration purposes
[nsresult, nsSegmentInfo] = ns_GetSegmentInfo(hfile, SegmentList(channel));
[nsresult, nsSegmentSourceInfo] = ns_GetSegmentSourceInfo(hfile, SegmentList(channel), 1);

% Load the first 100 waveforms on each selected channel
[nsresult, timestamps_wf, waveforms, sampleCount, unitIDs] = ns_GetSegmentData(hfile, SegmentList(channel), [1 : 1 : 100]);
clear timestamps_wf sampleCount;

clear nsSegmentInfo;
clear nsSegmentSourceInfo;

%
% Reduce the data set to only the first 150 seconds of data
%

NeuralLabels = strvcat(EntityInfo(NeuralList).EntityLabel);

for cChannel = 1 : 1 : length(channel),
    % Have to figure out which Neural entities correspond with the selected segment entities
    list = strmatch(EntityInfo(SegmentList(channel(cChannel))).EntityLabel, NeuralLabels, 'exact');
    
    % Retrieve the data
    [nsresult, NeuralInfo] = ns_GetNeuralInfo(hfile, NeuralList(list));
    [nsresult, NeuralData] = ns_GetNeuralData(hfile, NeuralList(list), 1, max([EntityInfo(NeuralList(list)).ItemCount]));

    if (totaltime == 150)
        ind = find(NeuralData <= 150);
    else
        ind = [1:1:size(NeuralData, 1) * size(NeuralData, 2)];
    end
    
    % Get the neural timestamps
    NeuralData = reshape(NeuralData, size(NeuralData, 1) * size(NeuralData, 2), 1);
    timestamps(cChannel) = {NeuralData(ind)};
    % Match the neural events with their unit ID
    temp = ones(length(ind) / length(list), 1) * [NeuralInfo(:).SourceUnitID];
    units(cChannel) = {temp(:)};
    % Remember how many neural events were found
    ItemCount(cChannel) = length(timestamps{cChannel});
    
    NeuralData(:) = [];
end
clear NeuralData SegmentLabels NeuralLabels NeuralInfo;
EntityInfo = rmfield(EntityInfo, 'ItemCount');

% Close data file. Should be done by the library but just in case. 
ns_CloseFile(hfile);

% Unload DLL
clear mexprog;

for cChannel = 1 : 1 : length(channel),    
    % Plot neural activity on that channel
    figure;
    subplot(4, 1, 1);
    % Plot the first 100 waveforms on that channel and apply different
    % colors for different unitIDs
    % Unclassified units
    ind = find(unitIDs(:, cChannel) == 0);
    if (~isempty(ind))
        plot(waveforms(:, ind, cChannel), 'Color', [.8 .8 .8]); % in grey
        hold on;
    end
    % Classified Units
    Colors = ['y' 'm' 'c' 'r' 'g'];
    
    for cUnit = 1 : 1 : 5,
        ind = find(unitIDs(:, cChannel) == cUnit);
        if (~isempty(ind))
            plot(waveforms(:, ind, cChannel), Colors(cUnit));
            hold on;
        end
    end
 
    axis tight;
    title(['Waveforms - ' EntityInfo(SegmentList(channel(cChannel))).EntityLabel]);
    ylabel('Voltage (uV)');
    
    subplot(4, 1, 2);
    % This creates lines for the unit 
    a = repmat(timestamps{cChannel}', 2, 1);          % X Values don't change so just two rows with the same values
    b = repmat([0.6; 1.4], 1, ItemCount(cChannel));   % Start and end points of each line
    % This plots the lines for each unit
    % Unclassified units
    ind = find(units{cChannel} == 0);
    hold on;
    if (~isempty(ind))
        plot(a(:, ind), b(:, ind), 'Color', [.8 .8 .8]); % in grey
    end
    
    % Classified Units
    for cUnit = 1 : 1 : 5,
        ind = find(unitIDs(:, cChannel) == cUnit);
        if (~isempty(ind))
            plot(a(:, ind), b(:, ind), Colors(cUnit));
        end
    end
    
    title(['Raster plot - ' EntityInfo(SegmentList(channel(cChannel))).EntityLabel]);
    xlabel('Time [s]');
    set(gca, 'ytick', []);
    xlim([0 round(totaltime)]);
    clear a b;
    
    % Use a sliding window 100 ms (delta t) long (equivalent to linear kernel)
    dt = 0.1; % ms
    for i = totaltime : -stepsize : 0,
        gaussdata(round(i / stepsize + 1)) = sum(exp(-(i - timestamps{cChannel}(find((timestamps{cChannel} > i - 4 * dt) & (timestamps{cChannel} < i + 4 * dt)))').^2 / (2 * dt^2))) / sqrt(2 * 3.1415) / dt;
    end
    
    % Throw away firing rates higher than 200 because they are not real and show up at the beginning
    gaussdata(find(gaussdata > 200)) = 0;
    
    subplot(4, 1, 3);
    plot(time, gaussdata(1 : length(time)));
    title(['Approximated firing rate (Gaussian) - ' EntityInfo(SegmentList(channel(cChannel))).EntityLabel]);
    xlabel('Time [s]');
    ylabel('Rate [Hz]');
    clear gaussdata;
    
    % Interspike Interval Distribution (1 ms resolution)
    if (length(timestamps{cChannel}) > 1)  % Have to make sure that we have at least two data points
        nbins = 500; % ms
        dtime = (timestamps{cChannel}(2 : end) - timestamps{cChannel}(1 : end - 1)) * 1000;
        
        subplot(4, 1, 4);
        bar(histc(dtime, 1 : 1 : nbins));
        axis tight;
        title(['Interspike Interval Histogram - ' EntityInfo(SegmentList(channel(cChannel))).EntityLabel]);
        xlabel('Time [ms]');
        ylabel('# Spikes');
        clear dtime;
    end    
end
