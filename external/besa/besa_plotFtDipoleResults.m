function status = besa_plotFtDipoleResults(cfg, timelock_data, ...
    dipolefit_data)

% BESA_PLOTFTDIPOLERESULTS visualizes the dipole fitting results generated
% by the call of the function ft_dipolefitting in FieldTrip.
% Please note that this function requires BESA Plot to be installed.
%
% Parameters:
%     [cfg]
%         A structure containing the settings for the plotting.
%           [cfg.besaplot] 
%               Full path to the program BesaPlot. Under Windows:
%               cfg.besaplot = '"C:\Program Files (x86)\BESA\Plot\BesaPlot.exe"';
%           [cfg.filebasename]
%               Three additional files are going to be generated for BESA
%               Plot. This is the basename for all three files.
%           [cfg.datapath]
%               The path to the folder where the additional files for BESA
%               Plot are going to be generated.
%           [cfg.badchannels]
%               The number of bad channels in the structure timelock_data.
%
% 
%     [timelock_data]
%         A structure resulting from the function ft_timelockanalysis.
% 
%     [dipolefit_data]
%         A structure resulting from the function ft_dipolefitting.
% 
% 
% Return:
%     [status] 
%         The status of the plotting process: 1 if the process was 
%         successful and less than 1 if not.
% 

% Copyright (C) 2016, BESA GmbH
%
% File name: besa_plotFtDipoleResults.m
%
% This file is part of MATLAB2BESA.
%
%    MATLAB2BESA is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    MATLAB2BESA is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with MATLAB2BESA. If not, see <http://www.gnu.org/licenses/>.
%
% Author: Todor Jordanov
% Created: 2016-07-06

status = 1;

%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variable %
%%%%%%%%%%%%%%%%%%%%%%%

BesaPlotPath = cfg.besaplot;
DataFile = [cfg.filebasename '.avr'];
ElpFile = [cfg.filebasename '.elp'];
ControlFile = [cfg.filebasename '.bpctrl'];
DataFilePath = fullfile(cfg.datapath, DataFile);
ControlPath = fullfile(cfg.datapath, ControlFile);
DatasetPath = [DataFilePath '.dataset'];
NumBadChannels = cfg.badchannels;
% Get the number of dipoles
NumDipoles = size(dipolefit_data.dip.pos, 1);

% The type of the channels to be stored.
% only EEG supported until now in BESA Plot
channel_type = 'EEG';

NumTimeSamples = size(timelock_data.avg, 2);
FirstTimeSample = timelock_data.time(1)*1000;
SamplingInterval = (timelock_data.time(2) - timelock_data.time(1))*1000;
NumChannels = size(timelock_data.avg, 1);
SortedCoordinates = zeros(NumChannels, 3);

data_matrix = timelock_data.avg;
time_samples = timelock_data.time;
channel_labels = timelock_data.label;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the sensor level data as an AVR-file. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isfield(timelock_data, 'grad') && isfield(timelock_data.grad, 'chanunit') && timelock_data.grad.chanunit{1} == 'T')
    
    data_scale_factor = 1e15; % fT
    
else
    
%     data_scale_factor = 1e6; % uV
    data_scale_factor = 1;

end

time_scale_factor = 1000; % ms

if(status == 1)
    
    status = besa_save2Avr(cfg.datapath, DataFile, data_matrix, ...
        time_samples, channel_labels, data_scale_factor, time_scale_factor);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rearrange channels in grad %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(status == 1)
    
    for iCh1 = 1:NumChannels

        CurrLabel1 = timelock_data.label{iCh1};

        for iCh2 = 1:NumChannels+NumBadChannels

            CurrLabel2 = timelock_data.grad.label{iCh2};

            if(strcmp(CurrLabel1, CurrLabel2))

                SortedCoordinates(iCh1, :) = timelock_data.grad.chanpos(iCh2, :);

            end

        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform to spherical coordinates %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(status == 1)
    
    SphericalCoords = zeros(NumChannels, 3);

    % Create a matrix for rotation about the z-axis.
    Angle1 = -90;
    % Angle1 = 0;
    rotate1 = [cosd(Angle1) -sind(Angle1) 0; sind(Angle1) cosd(Angle1) 0; 0 0 1];

    % Perform rotation.
    RotatedPositions = SortedCoordinates*rotate1;

    for iCh = 1:NumChannels

        % Get the current coordinates and normalize the radius to 1.0.
        CurrCartesianCoords = RotatedPositions(iCh, :) / ...
            sqrt(sum(RotatedPositions(iCh, :).^2));
        % Perform the transformation from Cartesian to spherical coordinates.
        [azimuth, elevation, r] = besa_transformCartesian2Spherical( ...
            CurrCartesianCoords(1), CurrCartesianCoords(2), ...
            CurrCartesianCoords(3));
        % Assign the transform values to the output matrix.
        SphericalCoords(iCh, 1) = azimuth;
        SphericalCoords(iCh, 2) = elevation;
        SphericalCoords(iCh, 3) = r;

    end

end

%%%%%%%%%%%%%%%%%%%
% Export elp-file %
%%%%%%%%%%%%%%%%%%%

if(status == 1)
    
    status = besa_save2Elp(cfg.datapath, ElpFile, ...
        SphericalCoords, channel_labels, channel_type);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the control file. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(status == 1)
    
    % Open file for writing.
    fid = fopen(ControlPath, 'w');

    % MATLAB reserves file identifiers 0, 1, and 2 for standard input,  
    % standard output (the screen), and standard error, respectively. When 
    % fopen successfully opens a file, it returns a file identifier greater 
    % than or equal to 3.
    if(fid >= 3)

        fprintf(fid, 'set color=black\n');
        fprintf(fid, 'set size=16\n');

        TimeSample = 1;
        for i = 1:NumDipoles

            % Get the current dipole location
            CurrDipLocation = dipolefit_data.dip.pos(i, :);
            % Shift 4 cm downwards
    %         CurrDipLocation(3) = CurrDipLocation(3) - 4.0;
            % Normalize
            NormedDipLoc = CurrDipLocation ./ sqrt(sum((CurrDipLocation).^2));
            NormedDipLoc(3) = NormedDipLoc(3) - 4.0/sqrt(sum((CurrDipLocation).^2));
            % Get current orientation
            CurrDipOri = dipolefit_data.dip.mom(3*i-2:3*i, TimeSample);
            % Normalize
            NormedDipOri = CurrDipOri ./ sqrt(sum((CurrDipOri).^2));
            fprintf(fid, 'set dipole number=%i xloc=%f yloc=%f zloc=%f xori=%f yori=%f zori=%f color=blue length=0.1 diameter=0.05\n', i, NormedDipLoc(2), NormedDipLoc(1), NormedDipLoc(3), NormedDipOri(2), NormedDipOri(1), NormedDipOri(3));

        end

        fprintf(fid, 'head width=35 height=40 viewpoint=top brain=on x=20.0 y=75.0\n');
        fprintf(fid, 'head width=35 height=40 viewpoint=back brain=on x=20.0 y=30.0\n');
        fprintf(fid, 'head width=35 height=40 viewpoint=left brain=on x=77.0 y=77.0\n');

        % Plot map
        % Average over time
        IndexOfTimeStart = find(timelock_data.time(:) == dipolefit_data.time(1));
        IndexOfTimeEnd = find(timelock_data.time(:) == dipolefit_data.time(end));
        avg_over_time = mean(timelock_data.avg(:, IndexOfTimeStart:IndexOfTimeEnd), 2);
        MapMax = round((max(avg_over_time).*data_scale_factor));
        MapMin = round((min(avg_over_time).*data_scale_factor));
        ContourStep = round((MapMax - MapMin) / 50.0);
        fprintf(fid, 'set dataset=1\n');
%         fprintf(fid, 'set mapchanneltype=Magnetometer\n');
        fprintf(fid, 'set mapuV=%i\n', ContourStep);
        fprintf(fid, 'set mapcolor negative=blue positive=red electrode=none\n');
        fprintf(fid, 'set mapsize=35\n');
        fprintf(fid, 'set maptype=amplitude\n');
        fprintf(fid, 'set maplambda=1E-5\n');
        fprintf(fid, 'set mapviewpoint=top\n');
        fprintf(fid, 'set maprange maximumuV=%i minimumuV=%i zerouV=0.0 contours=on\n', MapMax, MapMin);
        fprintf(fid, 'map latency=%.1f end=%.1f step=20 mean=on label=on x=77.0 y=30.0\n', dipolefit_data.time(1)*1000, dipolefit_data.time(end)*1000);
%         fprintf(fid, 'mapscale unit=1 side=left skip=5 x=55.0 y=30.0\n');

        fclose(fid);

    else

        status = 0;
        disp('Error! Invalid control file identifier.')

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the dataset file %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if(status == 1)
    
    status = besa_save2Dataset(DatasetPath, DataFile, ElpFile, ...
        ControlFile, NumChannels, NumTimeSamples, SamplingInterval, ...
        FirstTimeSample);

end

%%%%%%%%%%%%%%%%%%%%%%%
% Plot with BESA Plot %
%%%%%%%%%%%%%%%%%%%%%%%

if(status == 1)
    
    % Start BESA Plot in background.
    BesaPlotCommand = ['start /b cmd /c ' BesaPlotPath ' ' DataFilePath ' ' ControlPath];
    system(BesaPlotCommand, '-echo');

end
