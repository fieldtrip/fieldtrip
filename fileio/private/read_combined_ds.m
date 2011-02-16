function [dat] = read_combined_ds(dirname, hdr, begsample, endsample, chanindx)

% READ_COMBINED_DS
%
% Use as
%   hdr = read_combined_ds(dirname)
%   dat = read_combined_ds(dirname, hdr, begsample, endsample, chanindx)

% Copyright (C) 2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

needhdr = nargin==1 || isempty(hdr);
needdat = nargin>1;

supported = {
  'plexon_nex'
  'neuralynx_bin'
  'neuralynx_ncs'
  'fcdc_matbin'
  };

if needhdr

  ls = dir(dirname);                % get the list of filenames
  ls = ls(~cell2mat({ls.isdir}));   % remove the subdirs

  nfiles = length(ls);
  fname = cell(nfiles,1);
  ftype = cell(nfiles,1);
  sel   = false(nfiles,1);

  for i=1:nfiles
    fname{i} = fullfile(dirname, ls(i).name);
    ftype{i} = ft_filetype(fname{i});
    sel(i)   = any(strcmp(ftype{i}, supported));
    [p, f, x] = fileparts(fname{i});
    if filetype_check_extension(fname{i}, '.mat')
      % select only the *.bin and not the *.mat of each pair
      sel(i) = false;
    end
  end

  if ~any(sel)
    error('no supported files were found');
  end

  fname = fname(sel);
  ftype = ftype(sel);
  nfiles = length(fname);

  clear filehdr
  for i=1:nfiles
    filehdr(i) = ft_read_header(fname{i}, 'headerformat', ftype{i});
  end

  if any([filehdr.nChans]~=1)
    error('more than one channel per file not supported');
  else
    hdr.nChans = sum([filehdr.nChans]);
  end

  if length(unique([filehdr.label]))~=nfiles
    error('not all channels have a unique name');
  else
    hdr.label = [filehdr.label]';
  end

  if any(diff([filehdr.Fs]))
    error('different sampling frequenties per file not supported');
  else
    hdr.Fs = filehdr(1).Fs;
  end

  if any(diff([filehdr.nSamples]))
    error('different nSamples per file not supported');
  else
    hdr.nSamples = filehdr(1).nSamples;
  end

  if any(diff([filehdr.nSamplesPre]))
    error('different nSamplesPre per file not supported');
  else
    hdr.nSamplesPre = filehdr(1).nSamplesPre;
  end

  if any(diff([filehdr.nTrials]))
    error('different nTrials per file not supported');
  else
    hdr.nTrials = filehdr(1).nTrials;
  end

  if isfield(filehdr, 'TimeStampPerSample') && any(diff([filehdr.TimeStampPerSample]))
    error('different TimeStampPerSample per file not supported');
  else
    hdr.TimeStampPerSample = filehdr(1).TimeStampPerSample;
  end

  if isfield(filehdr, 'FirstTimeStamp') && any(diff([filehdr.FirstTimeStamp]))
    error('different FirstTimeStamp per file not supported');
  else
    hdr.FirstTimeStamp = filehdr(1).FirstTimeStamp;
  end

  % remember the original header details
  hdr.orig.header = filehdr;
  hdr.orig.fname  = fname;
  hdr.orig.ftype  = ftype;
end

if needdat
  filehdr = hdr.orig.header;
  fname   = hdr.orig.fname;
  ftype   = hdr.orig.ftype;

  if nargin<3
    begsample=1;
  end
  if nargin<4
    endsample = hdr.nSamples*hdr.nTrials;
  end
  if nargin<5
    chanindx = 1:hdr.nChans;
  end

  nsamples  = endsample-begsample+1;
  nchans    = length(chanindx);
  dat       = zeros(nchans, nsamples);

  for i=1:nchans
    thischan = chanindx(i);
    tmp = read_data(fname{thischan}, 'header', filehdr(thischan), 'dataformat', ftype{thischan}, 'begsample', begsample, 'endsample', endsample);
    dat(i,:) = tmp;
  end
end

if ~needdat
  % return the header
  dat = hdr;
end

