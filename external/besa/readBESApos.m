function sensor = readBESApos(filename, type)

% readBESApos read sensor information from a *.pos file. BESA stores
% magnetometer information in *.pmg, and gradiometer information in *.pos,
% however BESA Research doesn’t mind which extension is used, the
% distinction between gradiometers and magnetometers is based on the number
% of values on each line in the file.
% Brief format description:
% - one sensor per line
% - magnetometers: label (optional), six coordinates per line (location,
%   orientation)
% - gradiometers: label (optional), nine coordinates per line (location of
%   primary sensor, location of secondary sensor, orientation).
%
% Parameters:
%     [filename]
%         In the case that the current folder is not the folder containing 
%         the file it should be the full path including the name of the 
%         pos file else only the name of the file should be specified.
%
%     [type]
%         Optional parameter defining specific type of sensors that should
%         be read from file. 
%         Possible values: 'all' for all types of sensors
%                          'mag' for magnetometer sensors only
%                          'grad' for gradiometer sensors only
%         Default value is 'all'.
% 
% Return:
%     [sensor] 
%         The output is a structure containing the following fields:
%               NumSensors:     Defines number of sensors.
%               Label:          Label of each sensor.
%               Type:           Type of each sensor (either mag or grad).
%               Coordinates:    Coordinates of each sensor. Three values
%                               for magnetometers. Six values (location of
%                               primary sensor & location of secondary
%                               sensor) for gradiometers.
%               Orientations:   Orientations of each sensor.
% 
% Copyright (C) 2014, BESA GmbH
%
% File name: readBESApos.m
%
% Author: Robert Spangler
% Created: 2014-01-08

if isempty(findstr(filename,'.'))
  filename = [filename,'.pos'];
end

if exist('type', 'var')
    % Check if valid type defined
    if (strcmp(type, 'all') ~= 1 &&...
        strcmp(type, 'mag') ~= 1 &&...
        strcmp(type, 'grad') ~= 1)
            error('Invalid sensor type.');    
    end;
else
    % Set to default type
    type = 'all';
end;

fid = fopen(filename);

% MATLAB reserves file identifiers 0, 1, and 2 for standard input,  
% standard output (the screen), and standard error, respectively. When 
% fopen successfully opens a file, it returns a file identifier greater 
% than or equal to 3.
if(fid >= 3)
    % Number of read lines so far
    LineCnt = 0;
    % Number of sensors
    SensorCnt = 0;
    
    % Read line by line and store information according to number of values
    % in each line:
    % 6 values: magnetometer, values 1-3: (x/y/z) coordinates
    %                         values 4-6: orientations
    % 7 values: magnetometer, value 1: label;
    %                         values 1-3: (x/y/z) coordinates
    %                         values 4-6: orientations
    % 9 values: gradiometer, values 1-3: (x/y/z) coords of primary sensor
    %                        values 4-6: (x/y/z) coords of secondary sensor
    %                        values 7-9: orientations
    % 10 values: gradiometer, value 1: label
    %                         values 2-4: (x/y/z) coords of primary sensor
    %                         values 5-7: (x/y/z) coords of secondary sensor
    %                         values 8-10: orientations
    
    while(true)
        % Read line
        LineCnt = LineCnt + 1; 
        currline = fgetl(fid);
        if ~ischar(currline)
            break;
        end;
        
        % Separate number of values in current line
        tmp = regexp(currline, '\s+', 'split');
        % Remove empty cells
        tmp(cellfun(@isempty, tmp)) = [];
        % Number of entries in tmp cell array
        NumEntries = size(tmp, 2);
        if NumEntries == 6
            % Magnetometer
            currlabel = '';
            currtype = 'mag';
            currcoordinates = ...
                [str2double(tmp{1}) str2double(tmp{2}) str2double(tmp{3})];
            currorientations = ...
                [str2double(tmp{4}) str2double(tmp{5}) str2double(tmp{6})];
        elseif NumEntries == 7
            % Magnetometer
            currlabel = tmp{1};
            currtype = 'mag';
            currcoordinates = ...
                [str2double(tmp{2}) str2double(tmp{3}) str2double(tmp{4})];
            currorientations = ...
                [str2double(tmp{5}) str2double(tmp{6}) str2double(tmp{7})];
        elseif NumEntries == 9
            % Gradiometer
            currlabel = '';
            currtype = 'grad';
            currcoordinates = ...
                [str2double(tmp{1}) str2double(tmp{2}) str2double(tmp{3})...
                 str2double(tmp{4}) str2double(tmp{5}) str2double(tmp{6})];
            currorientations = ...
                [str2double(tmp{7}) str2double(tmp{8}) str2double(tmp{9})];
        elseif NumEntries == 10
            % Gradiometer
            currlabel = tmp{1};
            currtype = 'grad';
            currcoordinates = ...
                [str2double(tmp{2}) str2double(tmp{3}) str2double(tmp{4})...
                 str2double(tmp{5}) str2double(tmp{6}) str2double(tmp{7})];           
            currorientations = ...
                [str2double(tmp{8}) str2double(tmp{9}) str2double(tmp{10})];
        else
            % Number of values does not match number of entries for any of
            % the known formats.
            error('Inconsistent number of values in line.');
        end;
        
        % Store info in sensor struct if current sensor type matches
        % desired type of sensors ('all', ''mag, 'grad')        
        if strcmp(type, 'all') == 1
            addsensor = true;
        elseif strcmp(currtype, type) == 1
            addsensor = true;
        else
            % Neglect this sensor
            addsensor = false;
        end;
        if(addsensor)
            SensorCnt = SensorCnt + 1;
            sensor.Label{SensorCnt, 1} = currlabel;
            sensor.Type{SensorCnt, 1} = currtype;
            sensor.Coordinates{SensorCnt, :} = currcoordinates;
            sensor.Orientations(SensorCnt,:) = currorientations;   
        end;
    end;   

    % Close file
    fclose(fid);
    
    % store number of sensors in sensor struct
    sensor.NumSensors = SensorCnt;
else
    error('Unable to open file!');
end;

