function [vhdr] = read_brainvision_vhdr(filename)

% READ_BRAINVISION_VHDR reads the known items from the BrainVision EEG
% header file and returns them in a structure.
%
% Use as
%   vhdr = read_brainvision_vhdr(filename)
%
% See also READ_BRAINVISION_EEG, READ_BRAINVISION_VMRK

% Copyright (C) 2003-2019, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

vhdr.DataFile         = read_asa(filename, 'DataFile=', '%s');
vhdr.MarkerFile       = read_asa(filename, 'MarkerFile=', '%s');
vhdr.DataFormat       = read_asa(filename, 'DataFormat=', '%s');
vhdr.DataOrientation  = read_asa(filename, 'DataOrientation=', '%s');
vhdr.BinaryFormat     = read_asa(filename, 'BinaryFormat=', '%s');
vhdr.NumberOfChannels = read_asa(filename, 'NumberOfChannels=', '%d');
vhdr.SamplingInterval = read_asa(filename, 'SamplingInterval=', '%f'); % microseconds

if ~isempty(vhdr.NumberOfChannels)
  for i=1:vhdr.NumberOfChannels
    chan_str  = sprintf('Ch%d=', i);
    chan_info = read_asa(filename, chan_str, '%s');
    t = tokenize(chan_info, ',');
    vhdr.label{i} = t{1};
    vhdr.reference{i} = t{2};
    resolution = str2num(t{3}); % in microvolt
    if ~isempty(resolution)
      vhdr.resolution(i) = resolution;
    else
      ft_warning('unknown resolution (i.e. recording units) for channel %d in %s', i, filename);
      vhdr.resolution(i) = 1;
    end
  end
end

% compute the sampling rate in Hz
vhdr.Fs = 1e6/(vhdr.SamplingInterval);

% the number of samples is unkown to start with
vhdr.nSamples = Inf;

% confirm the names of the .vmrk and .eeg files as specified in the .vhdr file
[p, f, x] = fileparts(filename);
dataFile = fullfile(p, vhdr.DataFile);
[dum, dum, dataExt] = fileparts(vhdr.DataFile);

if ~isempty(vhdr.MarkerFile)
  markerFile = fullfile(p, vhdr.MarkerFile);
  [dum, dum, markerExt] = fileparts(vhdr.MarkerFile);
else
  % it is allowed to have a BrainVision dataset with only header and data, but no marker file
  markerFile = '';
  markerExt = '';
end

% in case they cannot be found, we can try whether they exist with the same name as the vhdr file
% this can happen if the researcher renamed the three files but not the header information
dataFileSameName = fullfile(p,[f dataExt]);
markerFileSameName = fullfile(p,[f markerExt]);

info = dir(dataFile);
if isempty(info)
  info = dir(filename);
  if ~isempty(info)
    ft_notice('Could not find "%s" as specified in the .vhdr file, using "%s" instead', dataFile, dataFileSameName);
    vhdr.DataFile = [f dataExt]; % without full path
    dataFile = dataFileSameName; % with full path
  else
    ft_error('cannot determine the location of the data file %s', dataFile);
  end
else
  if ~strcmp(dataFile, dataFileSameName)
    ft_notice('Could not find "%s" as specified in the .vhdr file, using "%s" instead', dataFile, dataFileSameName);
  end
end

if ~isempty(markerFile)
  info = dir(markerFile);
  if isempty(info)
    info = dir(markerFileSameName);
    if ~isempty(info)
      ft_notice('Could not find "%s" as specified in the .vhdr file, using "%s" instead', markerFile, markerFileSameName);
      vhdr.MarkerFile = [f markerExt]; % without full path
      markerFile = markerFileSameName; % with full path
    else
      ft_error('cannot determine the location of the marker file %s', markerFile);
    end
  else
    if ~strcmp(markerFile, markerFileSameName)
      ft_notice('Could not find "%s" as specified in the .vhdr file, using "%s" instead', markerFile, markerFileSameName);
    end
  end
end

% determine the number of samples by looking at the binary file
if strcmpi(vhdr.DataFormat, 'binary')
  % the data file is supposed to be located in the same directory as the header file
  % but that might be on another location than the present working directory
  
  info = dir(dataFile);
  switch lower(vhdr.BinaryFormat)
    case 'int_16'
      vhdr.nSamples = info.bytes./(vhdr.NumberOfChannels*2);
    case 'int_32'
      vhdr.nSamples = info.bytes./(vhdr.NumberOfChannels*4);
    case 'ieee_float_32'
      vhdr.nSamples = info.bytes./(vhdr.NumberOfChannels*4);
  end
  
elseif strcmpi(vhdr.DataFormat, 'ascii')
  vhdr.skipLines = 0;
  vhdr.skipColumns = 0;
  
  % Read ascii info from header (if available).
  dataPoints = read_asa(filename, 'DataPoints=', '%d');
  skipLines = read_asa(filename, 'SkipLines=', '%d');
  skipColumns = read_asa(filename, 'SkipColumns=', '%d');
  decimalSymbol = read_asa(filename, 'DecimalSymbol=', '%s'); % This is not used in reading dataset yet
  
  if ~isempty(dataPoints); vhdr.nSamples = dataPoints; end
  if ~isempty(skipLines); vhdr.skipLines = skipLines; end
  if ~isempty(skipColumns); vhdr.skipColumns = skipColumns; end
  if ~isempty(decimalSymbol); vhdr.decimalSymbol = decimalSymbol; end
  
  if isempty(dataPoints) && strcmpi(vhdr.DataOrientation, 'vectorized')
    % this is a very inefficient fileformat to read data from, it looks like this:
    % Fp1   -2.129 -2.404 -18.646 -15.319 -4.081 -14.702 -23.590 -8.650 -3.957
    % AF3   -24.023 -23.265 -30.677 -17.053 -24.889 -35.008 -21.444 -15.896 -12.050
    % F7    -10.553 -10.288 -19.467 -15.278 -21.123 -25.066 -14.363 -10.774 -15.396
    % F3    -28.696 -26.314 -35.005 -27.244 -31.401 -39.445 -30.411 -20.194 -16.488
    % FC1   -35.627 -29.906 -38.013 -33.426 -40.532 -49.079 -38.047 -26.693 -22.852
    % ...
    fid = fopen_or_error(dataFile, 'rt');
    tline = fgetl(fid); % read the complete first line
    fclose(fid);
    t = tokenize(tline, ' ', true); % cut the line into pieces
    vhdr.nSamples = length(t) - 1; % the first element is the channel label
  end
end

if isinf(vhdr.nSamples)
  ft_warning('cannot determine number of samples for this sub-fileformat');
end

% the number of trials is unkown, assume continuous data
vhdr.nTrials = 1;
vhdr.nSamplesPre = 0;

% ensure that the labels are in a column
vhdr.label = vhdr.label(:);
vhdr.reference = vhdr.reference(:);
vhdr.resolution = vhdr.resolution(:);

% read in impedance values
vhdr.impedances.channels = [];
vhdr.impedances.reference = [];
vhdr.impedances.ground = NaN;
vhdr.impedances.refChan = [];

fid = fopen_or_error(filename, 'rt');
while ~feof(fid)
  tline = fgetl(fid);
  if startsWith(tline, 'Impedance [')
    chanCounter = 0;
    refCounter = 0;
    impCounter = 0;
    while chanCounter<vhdr.NumberOfChannels && ~feof(fid)
      chan_info = fgetl(fid);
      if ~isempty(chan_info)
        impCounter = impCounter+1;
        [chanName,impedances] = strtok(chan_info,':');
        spaceList = strfind(chanName,' ');
        if ~isempty(spaceList)
          chanName = chanName(spaceList(end)+1:end);
        end
        if strfind(chanName,'REF_') == 1 %for situation where there is more than one reference
          refCounter = refCounter+1;
          vhdr.impedances.refChan(refCounter) = impCounter;
          if ~isempty(impedances)
            vhdr.impedances.reference(refCounter) = str2double(impedances(2:end));
          else
            vhdr.impedances.reference(refCounter) = NaN;
          end
        elseif strcmpi(chanName,'ref') %single reference
          refCounter = refCounter+1;
          vhdr.impedances.refChan(refCounter) = impCounter;
          if ~isempty(impedances)
            vhdr.impedances.reference(refCounter) = str2double(impedances(2:end));
          else
            vhdr.impedances.reference(refCounter) = NaN;
          end
        else
          chanCounter = chanCounter+1;
          if ~isempty(impedances)
            vhdr.impedances.channels(chanCounter,1) = str2double(impedances(2:end));
          else
            vhdr.impedances.channels(chanCounter,1) = NaN;
          end
        end
      end
    end
    if ~feof(fid)
      tline='';
      while ~feof(fid) && isempty(tline)
        tline = fgetl(fid);
      end
      if ~isempty(tline)
        if strcmp(tline(1:4),'Ref:')
          refCounter = refCounter+1;
          [chanName,impedances] = strtok(tline,':');
          if ~isempty(impedances)
            vhdr.impedances.reference(refCounter) = str2double(impedances(2:end));
          else
            vhdr.impedances.reference(refCounter) = NaN;
          end
        end
        if strcmpi(tline(1:4),'gnd:')
          [chanName,impedances] = strtok(tline,':');
          vhdr.impedances.ground = str2double(impedances(2:end));
        end
      end
    end
    if ~feof(fid)
      tline='';
      while ~feof(fid) && isempty(tline)
        tline = fgetl(fid);
      end
      if ~isempty(tline)
        if strcmpi(tline(1:4),'gnd:')
          [chanName,impedances] = strtok(tline,':');
          vhdr.impedances.ground = str2double(impedances(2:end));
        end
      end
    end
  end
end
fclose(fid);
