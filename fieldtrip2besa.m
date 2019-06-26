function fieldtrip2besa(filename, data, varargin)

% FIELDTRIP2BESA saves a FieldTrip data structures to a corresponding BESA file. This
% export function is based on documentation that was provided by Todor Jordanov of
% BESA.
%
% Use as
%   fieldtrip2besa(filename, data)
% with data as obtained from FT_PREPROCESSING to export single trial data as a
% set of .avr files.
%
% Use as
%   fieldtrip2besa(filename, elec)
% or
%   fieldtrip2besa(filename, grad)
% with an electrode structure as obtained from FT_READ_SENS to export channel
% positions to an .elp file.
%
% Additional key-value pairs can be specified according to
%   channel = cell-array, can be used to make subset and to reorder the channels
%
% See also FIELDTRIP2SPSS, FIELDTRIP2FIFF

% parse the optional input arguments
channel = ft_getopt(varargin, 'channel');

% this requires the "MATLAB to BESA Export functions"
% which are available from http://www.besa.de/downloads/matlab/
ft_hastoolbox('matlab2besa', 1);

datatype = ft_datatype(data);
switch datatype

  case 'raw'
    %% write raw data as *.avr
    assert(isempty(channel), 'channel selection and reordering is not yet supported');

    NumTrials = length(data.trial);
    channel_labels = data.label;
    data_scale_factor = 1.0;
    time_scale_factor = 1.0;
    for iTr = 1:NumTrials
      % Multiply by 1000 to get the time in milliseconds.
      time_samples = data.time{iTr}.*1000;
      % The file name where data should be written.
      file_name = sprintf('%s%03d.avr', filename, iTr);
      % Multiply by 1e15 to get the data in femtoTesla.
      data_matrix = data.trial{iTr}.*1e15; %FIXME
      % Save the data
      besa_save2Avr(custom_path, file_name, data_matrix, time_samples, channel_labels, data_scale_factor, time_scale_factor);
    end

  case 'timelock'
    %% write timelocked data as *.avr
    assert(isempty(channel), 'channel selection and reordering is not yet supported');

    if isfield(data, 'trial') && strcmp(getdimord(data, 'trial'), 'rpt_chan_time')
      [NumTrials, NumChans, NumSamp] = size(data.trial);

      % Multiply by 1000 to get the time in milliseconds.
      time_samples = data.time.*1000;
      channel_labels = data.label;
      data_scale_factor = 1.0;
      time_scale_factor = 1.0;
      for iTr = 1:NumTrials
        % The file name where data should be written.
        file_name = sprintf('%s%03d.avr', filename, iTr);
        % Multiply by 1e15 to get the data in femtoTesla.
        data_matrix = reshape(data.trial(iTr, :, :), [NumChans NumSamp]).*1e15; % FIXME
        % Save the data
        besa_save2Avr(custom_path, file_name, data_matrix, time_samples, channel_labels, data_scale_factor, time_scale_factor);
      end

    elseif isfield(data, 'avg') && strcmp(getdimord(data, 'avg'), 'chan_time')
      % Multiply by 1000 to get the time in milliseconds.
      time_samples = data.time.*1000;
      channel_labels = data.label;
      data_scale_factor = 1.0;
      time_scale_factor = 1.0;
      % The file name where data should be written.
      file_name = sprintf('%s.avr', filename);
      % Multiply by 1e15 to get the data in femtoTesla.
      data_matrix = data.avg.*1e15; % FIXME
      % Save the data
      besa_save2Avr(custom_path, file_name, data_matrix, time_samples, channel_labels, data_scale_factor, time_scale_factor);

    else
      ft_error('unsupported data structure');
    end

  case {'elec', 'grad'}
    %% write channel data to *.elp

    channel_labels = data.label;
    NumChannels = length(data.label);

    % Rearrange channels in grad
    SortedCoordinates = zeros(NumChannels, 3);
    NumBadChannels = 1; % A106
    for iCh1 = 1:NumChannels
      CurrLabel1 = channel{iCh1};
      for iCh2 = 1:NumChannels+NumBadChannels
        CurrLabel2 = data.label{iCh2};
        if(strcmp(CurrLabel1, CurrLabel2))
          SortedCoordinates(iCh1, :) = data.chanpos(iCh2, :);
        end
      end
    end

    % Transform to spherical coordinates
    SphericalCoords = zeros(NumChannels, 3);
    % Create a matrix for rotation about the z-axis.
    Angle1 = -90;
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

    % The type of the channels to be stored.
    switch datatype
      case 'grad'
        channel_type = 'MEG';
      case 'elec';
        channel_type = 'EEG';
    end

    % Export elp-file
    besa_save2Elp(custom_path, filename, SphericalCoords, channel_labels, channel_type);

  otherwise
    ft_error('unsupported data structure');

end % switch type
