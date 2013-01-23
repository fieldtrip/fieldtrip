function ExampleAnalog()
% function OutputAnalogVector = ExampleAnalog()

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
end

if (cAnalog == 0)
    disp('No analog entities available!');
end

if (cEvent == 0)
    disp('No event entities available!');
end

% prompt to get the number of points
max_count = input('How many data points (x-axis) do you want max? ');

% Have user pick a channels or channels for further analysis
% Show the user how many channels are available
disp(' ');
disp(['There are ' num2str(cAnalog) ' analog channels.']);
disp(' ');

yes = input('Show first all analog channels? y/n ', 's');
if (yes == 'y')
    for i = 1 : 1 : length(AnalogList)
        chan = AnalogList(i);
        count = min(max_count, EntityInfo(chan).ItemCount);
        % Get the fist data points of the waveform and show it
        [nsresult, ContinuousCount, wave] = ns_GetAnalogData(hfile, AnalogList(i), 1, count);
        figure;
        plot(wave);
    end
    return;
end

disp(['First analog entity: ' num2str(AnalogList(1)) ]);
channel = input('Which data channels would you like to display? (e.g. 1 or [1 2 3])');
if (EntityInfo(channel).EntityType == 2)  % Have to check that the selected channel actually exists
else
    disp('Channel is not of type analog');
    return
end
clear cNeural cSegment;

% Throw away entity infos we don't need to save memory
EntityInfo = rmfield(EntityInfo, 'EntityType');

%
% Load the waveform data and do the analysis
%


count = min(max_count, EntityInfo(channel).ItemCount);

% Get the fist data points of the waveform and show it
[nsresult, ContinuousCount, wave] = ns_GetAnalogData(hfile, channel, 1, count);
plot(wave);
return;
