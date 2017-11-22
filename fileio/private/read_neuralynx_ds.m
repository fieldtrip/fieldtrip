function [dat] = read_neuralynx_ds(dirname, hdr, begsample, endsample, chanindx)

% READ_NEURALYNX_DS reads multiple single-channel Neuralynx files that are
% all contained in a single directory. Each file is treated as a single
% channel of a combined multi-channel dataset.
%
% Use as
%   [hdr] = read_neuralynx_ds(dirname)
%   [dat] = read_neuralynx_ds(dirname, hdr, begsample, endsample, chanindx)
%
% A Neuralynx dataset consists of a directory containing separate files,
% one for each channel. All Neuralynx datafiles starts with a 16k header
% (in ascii format), followed by an arbitrary number of data records. The
% format of the data records depend on the type of data contained in the
% channel (e.g. continuous or spike data).
%
% To read the timestamps of spike waveforms (nse) or clustered spikes (nts),
% the header should contain the fields
%   hdr.FirstTimeStamp
%   hdr.TimeStampPerSample
% These can only be obtained from the corresponding simultaneous LFP
% and/or MUA recordings.
%
% See also READ_NEURALYNX_NCS, READ_NEURALYNX_NSE, READ_NEURALYNX_NTS

% Copyright (C) 2006-2007, Robert Oostenveld
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

needhdr = (nargin==1);
needdat = (nargin>=2);

if needhdr
  % get the list of filenames
  ls = dir(dirname);
  ls = ls(~cell2mat({ls.isdir}));
  fname = {};
  for i=1:length(ls)
    fname{i} = fullfile(dirname, ls(i).name);
  end

  ftype = zeros(length(fname), 1);
  for i=1:length(fname)
    if     ft_filetype(fname{i}, 'neuralynx_ncs')
      ftype(i) = 1;
    elseif ft_filetype(fname{i}, 'neuralynx_nse')
      ftype(i) = 2;
    elseif ft_filetype(fname{i}, 'neuralynx_nts')
      ftype(i) = 3;
    end
  end

  % only remember the filenames that are relevant
  fname = fname(ftype>0);
  ftype = ftype(ftype>0);
  ftype_ncs = find(ftype==1);
  ftype_nse = find(ftype==2);
  ftype_nts = find(ftype==3);

  if length(fname)==0
    ft_error('the dataset directory contains no supported files');
  end

  for i=1:length(fname)
    % this will only work if all files within a dataset return a similar header structure
    switch ftype(i)
      case 1
        orig(i) = read_neuralynx_ncs(fname{i}, 0, 0);
      case 2
        orig(i) = read_neuralynx_nse(fname{i}, 0, 0);
      case 3
        orig(i) = read_neuralynx_nts(fname{i}, 0, 0);
      otherwise
        ft_error('unsupported');
    end
  end

  % combine the information from the different files in a single header
  for i=1:length(orig)
    if isfield(orig(i).hdr, 'NLX_Base_Class_Name')
      label{i}           = orig(i).hdr.NLX_Base_Class_Name;
    else
      label{i}           = orig(i).hdr.AcqEntName;
    end
    if isfield(orig(i).hdr, 'SubSamplingInterleave')
      SamplingFrequency(i) = orig(i).hdr.SamplingFrequency / orig(i).hdr.SubSamplingInterleave;
    else
      SamplingFrequency(i) = orig(i).hdr.SamplingFrequency;
    end
    ADBitVolts(i)        = orig(i).hdr.ADBitVolts;
    FirstTimeStamp(i)    = orig(i).hdr.FirstTimeStamp;
    LastTimeStamp(i)     = orig(i).hdr.LastTimeStamp;
    NRecords(i)          = orig(i).NRecords;
    % Note that the last timestamp corresponds with the first sample of the last
    % record and not with the last sample in the file.
  end

  for i=1:length(orig)
    % timestamps are measured in units of approximately 1us
    % in case of 32556 Hz sampling rate, there are approximately 30.7 timestamps per sample
    % in case of 1000 Hz sampling rate,  there are approximately 1000 timestamps per sample
    % note that the last timestamp in the original header corresponds with the
    % first sample of the last record, and not with the last sample
    switch ftype(i)
      case 1
        % ensure that the last timestamp refers to the last sample
        recordsize = 512; % each record contains 512 samples
      case 2
        recordsize = 32; % each record contains 32 samples
      case 3
        ft_error('this has not been implemented yet');
      otherwise
        ft_error('unsupported');
    end
    % ensure that the last timestamp refers to the last sample
    TimeStampPerSample(i) = double(LastTimeStamp(i)-FirstTimeStamp(i))/((NRecords(i)-1)*recordsize);  % this should be in double precision, since it can be fractional
    LastTimeStamp(i)      = LastTimeStamp(i) + uint64((recordsize-1)*TimeStampPerSample(i));          % this should be in uint64 precision
  end % for length(orig)

  if any(SamplingFrequency~=SamplingFrequency(1))
    ft_warning('not all channels have the same sampling rate');
  end

  if ~isempty(ftype_ncs)
    if any(FirstTimeStamp(ftype_ncs)~=FirstTimeStamp(ftype_ncs(1)))
      % there seems to be a matlab bug (observed in Matlab75 on windows) which causes this uint64 comparison to fail if there are exactly 8 files
      % therefore check once more after converting them to double
      if any(double(FirstTimeStamp(ftype_ncs))~=double(FirstTimeStamp(ftype_ncs(1))))
        ft_error('not all continuous channels start at the same time');
      end
    end
    if any(LastTimeStamp(ftype_ncs)~=LastTimeStamp(ftype_ncs(1)))
      % there seems to be a matlab bug (observed in Matlab75 on windows) which causes this uint64 comparison to fail if there are exactly 8 files
      % therefore check once more after converting them to double
      if any(double(LastTimeStamp(ftype_ncs))~=double(LastTimeStamp(ftype_ncs(1))))
        ft_warning('not all continuous channels end at the same time');
      end
    end
    if any(NRecords(ftype_ncs)~=NRecords(ftype_ncs(1)))
      ft_warning('not all continuous channels have the same number of records');
    end
  end % if ftype_ncs

  % construct the header that applies to all channels combined
  hdr.nChans         = length(label);
  hdr.label          = label;
  hdr.filename       = fname;
  hdr.nTrials        = 1;                           % it is continuous
  hdr.Fs             = SamplingFrequency(1);
  hdr.nSamplesPre    = 0;                           % it is continuous

  if ~isempty(ftype_ncs)
    % these elements are only relevant for continously sampled channels
    hdr.nSamples           = NRecords(1) * 512;
    hdr.FirstTimeStamp     = FirstTimeStamp(1);
    hdr.LastTimeStamp      = LastTimeStamp(1);
    hdr.TimeStampPerSample = TimeStampPerSample(1);
  end

  % remember the original header details
  hdr.orig = orig(:);

  % return the header
  dat = hdr;

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % read the data of the selected channels (i.e. files)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  if nargin<5
    % select all channels
    chanindx = 1:length(hdr.label);
  end
  nchan   = length(chanindx);
  nsample = endsample-begsample+1;
  dat     = zeros(nchan, nsample);

  for i=1:nchan
    thischan = chanindx(i);
    thisfile = hdr.filename{thischan};
    switch ft_filetype(thisfile)
    case 'neuralynx_ncs'
      % determine the records that contain the sample numbers of the requested segment
      begrecord  = ceil(begsample/512);
      endrecord  = ceil(endsample/512);
      ncs = read_neuralynx_ncs(thisfile, begrecord, endrecord);
      % copy the selected samples into the output
      begsel = begsample - (begrecord-1)*512;
      endsel = endsample - (begrecord-1)*512;
      dat(i,:) = ncs.dat(begsel:endsel);

    case 'neuralynx_nse'
      % read all spike waveforms and timestamps
      nse = read_neuralynx_nse(thisfile);
      % convert the timestamps into samples
      fprintf('%d timstamps\n', length(nse.TimeStamp));
      sample = double(nse.TimeStamp-hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
      sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
      dat(i,sample) = dat(i,sample) + 1;

    case 'neuralynx_nts'
      % read all timestamps
      nts = read_neuralynx_nts(thisfile);
      % convert the timestamps into samples
      sample = double(nse.TimeStamp-hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
      sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
      dat(i,sample) = dat(i,sample) + 1;

    end % switch ft_filetype
  end % for nchan
end % reading data

