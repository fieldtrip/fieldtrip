function [dat] = read_brainvision_eeg(filename, hdr, begsample, endsample);

% READ_BRAINVISION_EEG reads raw data from an EEG file
% and returns it as a Nchans x Nsamples matrix
%
% Use as
%   dat = read_brainvision_eeg(filename, hdr, begsample, endsample)
% where the header should be first read using read_brainvision_vhdr
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_VMRK

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_brainvision_eeg.m,v $
% Revision 1.2  2009/03/13 07:13:28  roboos
% added feedback to vectorized ascii subformat
% fixed problem in vectorized ascii subformat when lines would start with the channel label
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.5  2008/04/09 10:10:58  roboos
% added support for int_32
% renamed nChans into the original form
%
% Revision 1.4  2007/06/13 08:08:19  roboos
% changed single & into &&
%
% Revision 1.3  2004/03/30 11:47:50  roberto
% dos->unix, fixed bug in multiplexed binary
%
% Revision 1.2  2004/03/30 08:22:19  roberto
% fixed bug due to renaming NumberOfChannels -> nChans
%
% Revision 1.1  2004/03/30 07:23:28  roberto
% added to CVS repository, copyrights added and implemented multiple file
% formats
%

if strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'multiplexed') && strcmpi(hdr.BinaryFormat, 'int_16')
  fid = fopen(filename, 'rb', 'ieee-le');
  fseek(fid, hdr.NumberOfChannels*2*(begsample-1), 'cof');
  [dat, siz] = fread(fid, [hdr.NumberOfChannels, (endsample-begsample+1)], 'int16');
  fclose(fid);
  % compute real microvolts using the calibration factor (resolution)
  res = sparse(diag(hdr.resolution));
  dat = res * dat;

elseif strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'multiplexed') && strcmpi(hdr.BinaryFormat, 'int_32')
  fid = fopen(filename, 'rb', 'ieee-le');
  fseek(fid, hdr.NumberOfChannels*4*(begsample-1), 'cof');
  [dat, siz] = fread(fid, [hdr.NumberOfChannels, (endsample-begsample+1)], 'int32');
  fclose(fid);
  % compute real microvolts using the calibration factor (resolution)
  res = sparse(diag(hdr.resolution));
  dat = res * dat;

elseif strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'multiplexed') && strcmpi(hdr.BinaryFormat, 'ieee_float_32')
  fid = fopen(filename, 'rb', 'ieee-le');
  fseek(fid, hdr.NumberOfChannels*4*(begsample-1), 'cof');
  [dat, siz] = fread(fid, [hdr.NumberOfChannels, (endsample-begsample+1)], 'float32');
  fclose(fid);

elseif strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'vectorized') && strcmpi(hdr.BinaryFormat, 'ieee_float_32')
  fid = fopen(filename, 'rb', 'ieee-le');
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
  fid = fopen(filename, 'rt');
  for line=1:(begsample-1)
    % read first lines and discard the data in them
    str = fgets(fid);
  end
  dat = zeros(endsample-begsample+1, hdr.NumberOfChannels);
  for line=1:(endsample-begsample+1)
    str = fgets(fid);			% read a single line with Nchan samples
    str(find(str==',')) = '.';		% replace comma with point
    dat(line,:) = str2num(str);
  end
  fclose(fid);
  % transpose the data
  dat = dat';

elseif strcmpi(hdr.DataFormat, 'ascii') && strcmpi(hdr.DataOrientation, 'vectorized')
  % this is a very inefficient fileformat to read data from, since it requires to
  % read in all the samples of each channel and then select only the samples of interest
  fid = fopen(filename, 'rt');
  dat = zeros(hdr.NumberOfChannels, endsample-begsample+1);
  for chan=1:hdr.NumberOfChannels
    % this is very slow, so better give some feedback to indicate that something is happening
    fprintf('reading channel %d from ascii file to get data from sample %d to %d\n', chan, begsample, endsample);

    str = fgets(fid);             % read all samples of a single channel
    str(find(str==',')) = '.';		% replace comma with point

    if ~isempty(regexp(str(1:10), '[a-zA-Z]', 'once'))
      % the line starts with letters, not numbers: probably it is a channel label
      % find the first number and remove the preceding part
      sel   = regexp(str(1:10), ' [-0-9]');   % find the first number, or actually the last space before the first number
      label = str(1:(sel));                   % this includes the space
      str   = str((sel+1):end);               % keep only the numbers
    end
    
    % convert the string into numbers and copy the desired samples over
    % into the data matrix
    tmp = str2num(str);
    dat(chan,:) = tmp(begsample:endsample);
  end
  fclose(fid);

else
  error('unsupported sub-fileformat');
end

