function sfh = readBESAsfh(filename, chantype)
% READBESASFH reads BESA coregistration information from SFH-file. 
%
% Parameters:
%     [filename]
%         In the case that the current folder is not the folder containing 
%         the file it should be the full path including the name of the 
%         sfh file else only the name of the file should be specified. 
% 
%     [chantype]
%         Optional parameter. Defines a specific type of sensors that will
%         be stored in sfh. All other sensors will be ignored.
%         Options: 'all' - default value
%                  'fiducial' - only fiducials will be taken into account
%                  'hpi' - only head position indicators will be taken into
%                  account
%                  'eeg' - only eeg sensors will be taken into account
%
% Return:
%     [sfh] 
%         A Matlab structure containing the coregistration information
%         stored in the SFH-file.
% 
% Copyright (C) 2013, BESA GmbH
%
% File name: readBESAsfh.m
%
% Author: Todor Jordanov
% Created: 2013-10-02

IsReadingOK = 1;
DataLevel = 0;

% Open the sfh-file for reading.
fp = fopen(filename, 'r');

% Check if the file was opened succesfully.
if(fp < 0)
    
    disp('Error: File does not exist or not enough permissions!')
    IsReadingOK = 0;
    
end

% Set defaults
if ~exist('chantype', 'var'), chantype = 'all'; end

% Get the first line of the file. It looks something like that:
% NrOfPoints: 34
if(IsReadingOK == 1)

    % Read the first line.
    FirstLine = fgetl(fp);
    % First of all check if the parameter NrOfPoints exists.
    if(~isempty(strfind(FirstLine, 'NrOfPoints')))

        tmp = regexp(FirstLine, 'NrOfPoints: ', 'split');
        NumPoints = str2double(tmp{2});
        DataLevel = 1;

    else

        disp('Error: Unknown file format!')
        IsReadingOK = 0;

    end

end

% Read the surface points labels and coordinates
if(IsReadingOK == 1 && DataLevel == 1)

    NumStoredPts = 0;    
    
    for p = 1:NumPoints

        CurrentLine = fgetl(fp);
        
        % Check if end of file.
        if(~ischar(CurrentLine))

            IsReadingOK = 0;
            disp('Error: File is corrupt or has not a valid format!')
            break;

        end
        
        % Split the line to 8 components delimited with white spaces.
        tmp = regexp(CurrentLine, '[ ]*', 'split');
        
        % In the MEG case there is a leading white space.
        if(isempty(tmp{1}))
            
            tmp = tmp(2:end);
            
        end
        
        if(size(tmp) ~= 8)

            IsReadingOK = 0;
            disp('Error: Wrong number of surface points parameters!')
            break;

        end
        
        % Get type of current point
        switch tmp{1}(1:3)
            case 'Fid'
                % Fiducial
                strType = 'fiducial';
            case 'Ele'
                % Electrode
                strType = 'eeg';                
            case 'Meg'
                % MEG sensor
                strType = 'meg';                
            case 'HPI'
                % Head position indicator
                strType = 'hpi';
            otherwise
                strType = 'unknown';
        end

        % Check if current point needs to be stored.
        bStoreInfo = 0;
        if(strcmp(chantype, 'all') == 1)
            bStoreInfo = 1;
        elseif(strcmp(chantype, strType) == 1)
            bStoreInfo = 1;
        end

        % Store info on current point in corresponding cell array
        if(bStoreInfo == 1)
            NumStoredPts = NumStoredPts + 1;
            sfh.SurfacePointsLabels{NumStoredPts} = tmp{1};
            sfh.SurfacePointsCoordinates(NumStoredPts, :) = ...
                str2double(tmp(2:end));
            sfh.SurfacePointsTypes{NumStoredPts} = strType;
        end

    end
    
    % Set correct number of entries
    sfh.NrOfPoints = NumStoredPts;
    
    if(IsReadingOK == 1)
        
        DataLevel = 2;
        
    end

end

% Read the transformation into BV-coordinates.
if(IsReadingOK == 1 && DataLevel == 2)

    CurrentLine = fgetl(fp);
    % Check if end of file.
    if(~ischar(CurrentLine))

        disp('Warning: The file contains no coregistration information!')

    else
        
        % If the line is not empty it should be something like this:
        % # Trans-data(in BV-coords): 3 translation, 3 rotation (in grad), 3 scale
        % There is no important information in that line, hence read the
        % next one. It shoul contain nine double values.
        CurrentLine = fgetl(fp);
        %tmp = str2num(CurrentLine);
        tmp = sscanf(CurrentLine, '%f', [9,1]);
        if(length(tmp) < 9)
            
            disp('Error: Not enough transformation parameters!')
            IsReadingOK = 0;
            
        else
            
            sfh.Transformations.Translation = tmp(1:3);
            sfh.Transformations.Rotation = tmp(4:6);
            sfh.Transformations.Scale = tmp(7:9);
            DataLevel = 3;
        
        end
        
    end

end

% Read the fiducials information.
if(IsReadingOK == 1 && DataLevel == 3)
    
    % This line should contain only the following:
    % Fiducials:
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'Fiducials')))
        
        % The next three lines should contain the fiducial information.
        sfh.Fiducials = zeros(3, 3);
        for i = 1:3
            
            CurrentLine = fgetl(fp);
            sfh.Fiducials(i, :) = sscanf(CurrentLine, '%f', [3,1]);
            
        end
        
    else
        
        disp('Error: Missing fiducial information!')
        IsReadingOK = 0;
        
    end
    
    if(IsReadingOK == 1)
        
        DataLevel = 4;
        
    end
   
end

% Read the midpoint information.
if(IsReadingOK == 1 && DataLevel == 4)
    
    % This line should contain only the following:
    % Midpoint (in BV-coords):
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'Midpoint')))
        
        sfh.Midpoint = zeros(1, 3);
        CurrentLine = fgetl(fp);
        sfh.Midpoint(1, :) = sscanf(CurrentLine, '%f', [3,1]);
        
    else
        
        disp('Error: Missing midpoint information!')
        IsReadingOK = 0;
       
    end
    
    if(IsReadingOK == 1)
        
        DataLevel = 5;
        
    end
    
end

% Read the information about volume MRI data files.
if(IsReadingOK == 1 && DataLevel == 5)
    
    % Get volume data path.
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'Volume')))
        
        sfh.VolumeDataPath = sscanf(CurrentLine, 'Volume: %s');
        
    else
        
        disp('Error: No volume data file available!')
        IsReadingOK = 0;
        
    end
    
    if(IsReadingOK == 1)
        
        DataLevel = 6;
        
    end
    
end

% Read the path of surface MRI data file.
if(IsReadingOK == 1 && DataLevel == 6)
    
    % Get surface data path.
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'Surface')))
        
        sfh.SurfaceDataPath = sscanf(CurrentLine, 'Surface: %s');
        
    else
        
        disp('Error: No Surface data file available!')
        IsReadingOK = 0;
        
    end
    
    if(IsReadingOK == 1)
        
        DataLevel = 7;
        
    end
    
end

% Read Talairach transformation information.
if(IsReadingOK == 1 && DataLevel == 7)
    
    % The next 8 lines should contain the coordinates of the following points:
    % AC PC AP PP SP IP RP LP
    sfh.Talairach.Labels = {'AC' 'PC' 'AP' 'PP' 'SP' 'IP' 'RP' 'LP'};
    
    for j = 1:length(sfh.Talairach.Labels)
        
        CurrentLine = fgetl(fp);
        
        if(~isempty(strfind(CurrentLine, sfh.Talairach.Labels{j})))
        
            sfh.Talairach.Coords(j, :) = ...
                sscanf(CurrentLine, [sfh.Talairach.Labels{j} ': %f%f%f']);

        else

            disp('Error: Invalid Talairach point %i!', j)
            IsReadingOK = 0;

        end
        
    end
    
    if(IsReadingOK == 1)
        
        DataLevel = 8;
        
    end
    
end

% Read the path of Talairach volume data file.
if(IsReadingOK == 1 && DataLevel == 8)
    
    % Get volume data path.
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'TalVolume')))
        
        sfh.Talairach.TalVolumePath = sscanf(CurrentLine, 'TalVolume: %s');
        
    else
        
        disp('Error: No Talairach volume data file available!')
        IsReadingOK = 0;
        
    end
    
    if(IsReadingOK == 1)
        
        DataLevel = 9;
        
    end
    
end

% Read the path of Talairach surface data file.
if(IsReadingOK == 1 && DataLevel == 9)
    
    % Get surface data path.
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'TalSurface')))
        
        sfh.Talairach.TalSurfacePath = sscanf(CurrentLine, 'TalSurface: %s');
        
    else
        
        disp('Error: No Talairach surface data file available!')
        IsReadingOK = 0;
        
    end
    
    if(IsReadingOK == 1)
        
        DataLevel = 10;
        
    end
    
end

% Read the path of Talairach brain surface data file.
if(IsReadingOK == 1 && DataLevel == 10)
    
    % Get brain surface data path.
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'TalBrainSurface')))
        
        sfh.Talairach.TalBrainSurfacePath = sscanf(CurrentLine, 'TalBrainSurface: %s');
        
    else
        
        disp('Error: No Talairach brain surface data file available!')
        IsReadingOK = 0;
        
    end
    
    if(IsReadingOK == 1)
        
        DataLevel = 11;
        
    end
    
end

% Read the path to FEM lead field file.
if(IsReadingOK == 1 && DataLevel == 11)
    
    % Get FEM lead field file path.
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'FEMtable')))
        
        sfh.FEM.FEMtable = sscanf(CurrentLine, 'FEMtable: %s');
        DataLevel = 12;
        
    else
        
        disp('Warning: No FEM lead field file available!')
        
    end
    
end

% Read the path to FEM grid file.
if(IsReadingOK == 1 && DataLevel == 12)
    
    % Get FEM grid file path.
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'FEMgrid')))
        
        sfh.FEM.FEMgrid = sscanf(CurrentLine, 'FEMgrid: %s');
        DataLevel = 13;
        
    else
        
        disp('Warning: No FEM grid file available!')
        
    end
    
end

% Read the path to FEM surface file.
if(IsReadingOK == 1 && DataLevel == 13)
    
    % Get FEM surface file path.
    CurrentLine = fgetl(fp);
    
    if(~isempty(strfind(CurrentLine, 'FEMsurface')))
        
        sfh.FEM.FEMsurface = sscanf(CurrentLine, 'FEMsurface: %s');
        
    else
        
        disp('Warning: No FEM surface file available!')
        
    end
    
end

if(IsReadingOK == 1)
    
    fclose(fp);

end