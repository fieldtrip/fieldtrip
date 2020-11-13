function [dat] = read_brainvision_eeg(filename, hdr, begsample, endsample, chanindx)

% READ_BRAINVISION_EEG reads raw data from an EEG file
% and returns it as a Nchans x Nsamples matrix
%
% Use as
%   dat = read_brainvision_eeg(filename, hdr, begsample, endsample)
% where the header should be first read using read_brainvision_vhdr
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_VMRK

% Copyright (C) 2003-2011, Robert Oostenveld
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

if nargin<5
  % read all channels
  chanindx = [];
end

if isequal(chanindx(:)', 1:hdr.NumberOfChannels);
  % read all channels
  chanindx = [];
end

% FIXME it would be nice to also implement the efficient reading of the
% selected channels for the other file formats but at the moment only the
% implementation of the binary multiplexed formats is smart enough.

if strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'multiplexed') && any(strcmpi(hdr.BinaryFormat, {'int_16', 'int_32', 'ieee_float_32'}))
  
  switch lower(hdr.BinaryFormat)
    case 'int_16'
      sampletype = 'int16';
      samplesize = 2;
    case 'int_32'
      sampletype = 'int32';
      samplesize = 4;
    case 'ieee_float_32'
      sampletype = 'float32';
      samplesize = 4;
  end % case
  
  fid = fopen_or_error(filename, 'rb', 'ieee-le');
  
  if isempty(chanindx)
    % read all the channels
    fseek(fid, hdr.NumberOfChannels*samplesize*(begsample-1), 'cof');
    dat = fread(fid, [hdr.NumberOfChannels, (endsample-begsample+1)], sampletype);
    % compute real microvolts using the calibration factor (resolution)
    % calib = diag(hdr.resolution);
    % % using a sparse multiplication speeds it up
    % dat = full(sparse(calib) * dat);
    calib = reshape(hdr.resolution,[],1);
    for k = 1:size(dat,2)
      dat(:,k) = calib.*dat(:,k);
    end
    
  else
    % read only the selected channels
    dat = zeros(length(chanindx), endsample-begsample+1);
    for chan = length(chanindx):-1:1
      fseek(fid, hdr.NumberOfChannels*samplesize*(begsample-1) + (chanindx(chan)-1)*samplesize, 'bof');
      dat(chan,:) = fread(fid, [1, (endsample-begsample+1)], sampletype, (hdr.NumberOfChannels-1)*samplesize);
    end
    % compute real microvolts using the calibration factor (resolution)
    % calib = diag(hdr.resolution(chanindx));
    % % using a sparse multiplication speeds it up
    % dat = full(sparse(calib) * dat);
    calib = reshape(hdr.resolution(chanindx),[],1);
    for k = 1:size(dat,2)
      dat(:,k) = calib.*dat(:,k);
    end
    
    % don't do the channel selection again at the end of the function
    chanindx = [];
  end
  
  fclose(fid);
  
elseif strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'vectorized') && strcmpi(hdr.BinaryFormat, 'ieee_float_32')
  fid = fopen_or_error(filename, 'rb', 'ieee-le');
  fseek(fid, 0, 'eof');
  hdr.nSamples = ftell(fid)/(4*hdr.NumberOfChannels);
  fseek(fid, 0, 'bof');
  numsamples = (endsample-begsample+1);
  for chan=1:hdr.NumberOfChannels
    fseek(fid, (begsample-1)*4, 'cof');                 % skip the first N samples
    [tmp, siz] = fread(fid, numsamples, 'float32');     % read these samples
    fseek(fid, (hdr.nSamples-endsample)*4, 'cof');      % skip the last M samples
    dat(chan,:) = tmp(:)';
  end
  fclose(fid);
  
elseif strcmpi(hdr.DataFormat, 'ascii') && strcmpi(hdr.DataOrientation, 'multiplexed')
  fid = fopen_or_error(filename, 'rt');
  
  % skip lines if hdr.skipLines is set and not zero
  if isfield(hdr,'skipLines') && hdr.skipLines > 0
    for line=1:hdr.skipLines
      str = fgets(fid);
    end
  end
  
  for line=1:(begsample-1)
    % read first lines and discard the data in them
    str = fgets(fid);
  end
  dat = zeros(endsample-begsample+1, hdr.NumberOfChannels);
  for line=1:(endsample-begsample+1)
    str = fgets(fid);         % read a single line with Nchan samples
    str(str==',') = '.';      % replace comma with point
    dat(line,:) = str2num(str);
  end
  fclose(fid);
  % transpose the data
  dat = dat';
  
elseif strcmpi(hdr.DataFormat, 'ascii') && strcmpi(hdr.DataOrientation, 'vectorized')
  % this is a very inefficient fileformat to read data from, since it requires to
  % read in all the samples of each channel and then select only the samples of interest
  fid = fopen_or_error(filename, 'rt');
  dat = zeros(hdr.NumberOfChannels, endsample-begsample+1);
  skipColumns = 0;
  for chan=1:hdr.NumberOfChannels
    % this is very slow, so better give some feedback to indicate that something is happening
    fprintf('reading channel %d from ascii file to get data from sample %d to %d\n', chan, begsample, endsample);
    
    % check whether columns have to be skipped
    if isfield(hdr,'skipColumns'); skipColumns = hdr.skipColumns; end
    
    str = fgets(fid);             % read all samples of a single channel
    str(str==',') = '.';          % replace comma with point
    
    if ~isempty(regexp(str(1:10), '[a-zA-Z]', 'once'))
      % the line starts with letters, not numbers: probably it is a channel label
      % find the first number and remove the preceding part
      sel   = regexp(str(1:10), ' [-0-9]');   % find the first number, or actually the last space before the first number
      label = str(1:(sel));                   % this includes the space
      str   = str((sel+1):end);               % keep only the numbers
      
      % as heading columns are already removed, set skipColumns to zero
      skipColumns = 0;
    end
    
    % convert the string into numbers and copy the desired samples over
    % into the data matrix
    tmp = str2num(str);
    dat(chan,:) = tmp(skipColumns+begsample:skipColumns+endsample);
  end
  fclose(fid);
  
else
  ft_error('unsupported sub-fileformat');
end

if ~isempty(chanindx)
  % for the the multiplexed binary formats the channel selection was
  % already done in the code above, for the other formats the selection
  % still has to be done here
  dat = dat(chanindx,:);
end
