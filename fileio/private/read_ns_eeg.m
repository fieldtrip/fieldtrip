function [eeg] = read_ns_eeg(filename, epoch)

% READ_NS_EEG read a NeuroScan 3.x or 4.x EEG File
%
% [eeg] = read_ns_eeg(filename, epoch)
%
%   filename     input Neuroscan .eeg file (version 3.x)
%   epoch        which epoch to read (default is all)
% 
% The output data structure eeg has the fields:
%   eeg.data(..)    - epoch signal in uV (size: Nepoch x Nchan x Npnt)
% and
%   eeg.label       - electrode labels
%   eeg.nchan       - number of channels
%   eeg.npnt        - number of samplepoints in ERP waveform
%   eeg.time        - time for each sample
%   eeg.rate        - sample rate (Hz)
%   eeg.xmin        - prestimulus epoch start (e.g., -100 msec)
%   eeg.xmax        - poststimulus epoch end (e.g., 900 msec)
%   eeg.nsweeps     - number of accepted trials/sweeps

% Copyright (C) 2003, Robert Oostenveld
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

% read the neuroscan header 
eeg = read_ns_hdr(filename);

% clear the variance part which is empty anyway
eeg = rmfield(eeg, 'variance');

% create a time axis
eeg.time = linspace(eeg.xmin, eeg.xmax, eeg.npnt);

% open the file and seek towards the place where the raw data is
fid = fopen(filename,'r','ieee-le');
if fid<0
  error(['cannot open ', filename]);
end

% the default is to read all epochs
if nargin<2
  epoch = 1:eeg.nsweeps;
end

% determine whether it is 16 or 32 bit data
fseek(fid, 0, 'eof');
header_size = 900 + 75*eeg.nchan;
file_size   = ftell(fid);
sample_size = (file_size-header_size)/(eeg.nchan*eeg.npnt*eeg.nsweeps);
% note that the V4 format can have some extra information at the end of
% the file, causing the sample size to be slightly larger than 2 or 4
if floor(sample_size)==2
  epoch_size = eeg.nchan*eeg.npnt*2 + 13;
  datatype ='int16';
elseif floor(sample_size)==4
  datatype ='int32';
  epoch_size = eeg.nchan*eeg.npnt*4 + 13;
elseif floor(sample_size)>4
  % although this is not to be expected, it might be due to a very large "footer" compared to the data size
  % Olga Sysoeva reported that this extention to the datatype detection fixed it for her (see http://bugzilla.fcdonders.nl/show_bug.cgi?id=547)
  datatype ='int32';
  epoch_size = eeg.nchan*eeg.npnt*4 + 13;
end

% create empty storage for the data
data = zeros(length(epoch), eeg.nchan, eeg.npnt);
sweep = struct();

for i=1:length(epoch)
  fseek(fid, 900, 'bof');               % skip general header
  fseek(fid, 75*eeg.nchan, 'cof');          % skip channel headers
  status = fseek(fid, (epoch(i)-1)*epoch_size, 'cof');  % skip first epochs
  if status~=0
    error('seek error while reading epoch data');
  end

  % fprintf('reading epoch %d at offset %d\n', epoch(i), ftell(fid));
    
  % read sweep header   
  sweep(i).accept   = fread(fid, 1, 'uchar');
  sweep(i).type     = fread(fid, 1, 'ushort');
  sweep(i).correct  = fread(fid, 1, 'ushort');
  sweep(i).rt       = fread(fid, 1, 'float32');
  sweep(i).response = fread(fid, 1, 'ushort');
  sweep(i).reserved = fread(fid, 1, 'ushort');

  % read raw signal
  raw = fread(fid, eeg.nchan*eeg.npnt, datatype);
  if length(raw)~=eeg.nchan*eeg.npnt
    error('fread error while reading epoch data');
  end
  data(i,:,:) = reshape(raw, [eeg.nchan, eeg.npnt]);

end

% convert raw signal to uV
for chan=1:eeg.nchan
  data(:,chan,:) = (data(:,chan,:) - eeg.baseline(chan)) * eeg.factor(chan);
end

% store the epoch information and data in the output structure
eeg.data = squeeze(data);
eeg.sweep = squeeze(sweep');

fclose(fid);

