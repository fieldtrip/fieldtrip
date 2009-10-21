function [hdr] = read_brainvision_vhdr(filename);

% READ_BRAINVISION_VHDR reads the known items from the BrainVision EEG
% header file and returns them in a structure
%
% Use as
%   hdr = read_brainvision_vhdr(filename)
%
% See also READ_BRAINVISION_EEG, READ_BRAINVISION_VMRK

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_brainvision_vhdr.m,v $
% Revision 1.2  2009/03/12 17:09:31  roboos
% improved the warning in case number of samples cannot be determined (applies to ascii *.seg data)
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.7  2008/11/14 07:42:13  roboos
% use general tokenize function instead of local copy, removed tokenize as subfunction
%
% Revision 1.6  2008/07/10 07:15:30  roboos
% also determine data file size if on another directory, thanks to Paul
%
% Revision 1.5  2008/04/16 07:51:11  roboos
% added fixme comment
%
% Revision 1.4  2008/04/09 10:08:28  roboos
% added detection of number of samples, based on filesize (only if binary)
%
% Revision 1.3  2008/04/09 10:06:50  roboos
% renamed nChans into NumberOfChannels to stay as close as possible to the ascii header
%
% Revision 1.2  2004/03/30 11:23:50  roberto
% fixed bug in cutting channel info into pieces, added local tokenize function
%
% Revision 1.1  2004/03/30 07:23:28  roberto
% added to CVS repository, copyrights added and implemented multiple file
% formats
%

hdr.DataFile         = read_asa(filename, 'DataFile=', '%s');
hdr.MarkerFile       = read_asa(filename, 'MarkerFile=', '%s');
hdr.DataFormat       = read_asa(filename, 'DataFormat=', '%s');
hdr.DataOrientation  = read_asa(filename, 'DataOrientation=', '%s');
hdr.BinaryFormat     = read_asa(filename, 'BinaryFormat=', '%s');
hdr.NumberOfChannels = read_asa(filename, 'NumberOfChannels=', '%d');
hdr.SamplingInterval = read_asa(filename, 'SamplingInterval=', '%f');	% microseconds

if ~isempty(hdr.NumberOfChannels)
  for i=1:hdr.NumberOfChannels
    chan_str  = sprintf('Ch%d=', i);
    chan_info = read_asa(filename, chan_str, '%s');
    t = tokenize(chan_info, ',');
    hdr.label{i} = t{1};
    hdr.reference{i} = t{2};
    resolution = str2num(t{3});          % in microvolt
    if ~isempty(resolution)
      hdr.resolution(i) = resolution;
    else
      hdr.resolution(i) = nan;
    end
  end
end

% compute the sampling rate in Hz
hdr.Fs = 1e6/(hdr.SamplingInterval);

% the number of samples is unkown to start with
hdr.nSamples = Inf;

% determine the number of samples by looking at the binary file
if strcmp(hdr.DataFormat, 'BINARY')
  % the data file is supposed to be located in the same directory as the header file
  % but that might be on another location than the present working directory
  [p, f, x] = fileparts(filename);
  datafile = fullfile(p, hdr.DataFile);
  info = dir(datafile);
  if isempty(info)
    error('cannot determine the location of the data file %s', hdr.DataFile);
  end
  switch lower(hdr.BinaryFormat)
    case 'int_16';
      hdr.nSamples = info.bytes./(hdr.NumberOfChannels*2);
    case 'int_32';
      hdr.nSamples = info.bytes./(hdr.NumberOfChannels*4);
    case 'ieee_float_32';
      hdr.nSamples = info.bytes./(hdr.NumberOfChannels*4);
  end
end

if isinf(hdr.nSamples)
  warning('cannot determine number of samples for this sub-fileformat');
end

% the number of trials is unkown, assume continuous data
hdr.nTrials     = 1;
hdr.nSamplesPre = 0;

% ensure that the labels are in a column
hdr.label      = hdr.label(:);
hdr.reference  = hdr.reference(:);
hdr.resolution = hdr.resolution(:);
