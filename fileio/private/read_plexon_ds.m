function [dat] = read_plexon_ds(dirname, hdr, begsample, endsample, chanindx)

% READ_PLEXON_DS reads multiple single-channel Plexon files that are
% all contained in a single directory. Each file is treated as a single
% channel of a combined multi-channel dataset.
%
% Use as
%   hdr = read_plexon_ds(dirname)
%   dat = read_plexon_ds(dirname,  hdr, begsample, endsample, chanindx)
%
% See also READ_PLEXON_NEX, READ_PLEXON_PLX, READ_PLEXON_DDT

% Copyright (C) 2007, Robert Oostenveld
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
    % determine the filetype only based on the extension, not using the filetype function because that is much slower
    if     filetype_check_extension(fname{i}, '.nex')
      ftype(i) = 1;
    elseif filetype_check_extension(fname{i}, '.plx')
      ftype(i) = 2;
    elseif filetype_check_extension(fname{i}, '.ddt')
      ftype(i) = 3;
    end
  end

  % only remember the filenames that are relevant
  fname = fname(ftype>0);
  ftype = ftype(ftype>0);

  if isempty(fname)
    ft_error('the directory contains no supported files');
  elseif any(ftype~=1)
    ft_error('only nex files are supported in a plexon dataset directory');
  end

  for i=1:length(fname)
    % this will only work if all files within a dataset return a similar header structure
    switch ftype(i)
      case 1
        orig(i) = read_plexon_nex(fname{i});
      case 'plexon_plx'
        ft_error('plx files are not supported in plexon dataset directory');
      case 'plexon_ddt'
        ft_error('ddt files are not supported in plexon dataset directory');
      otherwise
        ft_error('unsupported file in plexon dataset directory');
    end
  end

  for i=1:length(fname)
    if length(orig(i).VarHeader)>1
      ft_error('multiple channels in a single NEX file not supported');
    else
      % combine the information from the different files in a single header
      label{i}        = deblank(orig(i).VarHeader.Name);
      Type(i)         = orig(i).VarHeader.Type;
      WFrequency(i)   = orig(i).VarHeader.WFrequency;  % of the waveform
      ADBitVolts(i)   = orig(i).VarHeader.ADtoMV;
      NPointsWave(i)  = orig(i).VarHeader.NPointsWave;
      Beg(i)          = orig(i).FileHeader.Beg;
      End(i)          = orig(i).FileHeader.End;
      Frequency(i)    = orig(i).FileHeader.Frequency;  % of the timestamps
    end
  end

  if any(Type~=5)
    ft_error('not all channels contain continuous data');
  end

  if any(WFrequency~=WFrequency(1))
    ft_warning('not all channels have the same sampling rate');
  end

  if any(Frequency~=Frequency(1))
    ft_warning('not all channels have the same timestamp rate');
  end

  if any(Beg~=Beg(1))
    ft_warning('not all channels start at the same time');
  end

  if any(End~=End(1))
    ft_warning('not all channels end at the same time');
  end

  if any(NPointsWave~=NPointsWave(1))
    ft_warning('not all channels have the same number of samples');
  end

  % construct the header that applies to all channels combined
  hdr.label               = label;
  hdr.nChans              = length(label);
  hdr.nSamples            = NPointsWave(1);
  hdr.nSamplesPre         = 0;                           % it is continuous
  hdr.nTrials             = 1;                           % it is continuous
  hdr.Fs                  = WFrequency(1);
  hdr.FirstTimeStamp      = Beg(1);
  hdr.LastTimeStamp       = End(1);                      % FIXME this is often not correct
  hdr.TimeStampPerSample  = Frequency(1)/WFrequency(1);
  hdr.filename            = fname;
  
  % remember the filename and filetype to speed up subsequent calls for reading the data
  hdr.filetype            = cell(size(fname));
  hdr.filetype(ftype==1)  = {'plexon_nex'};
  hdr.filetype(ftype==2)  = {'plexon_plx'};
  hdr.filetype(ftype==3)  = {'plexon_ddt'};
  
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
    thistype = hdr.filetype{thischan};
    
    switch thistype
      case 'plexon_nex'
        nex = read_plexon_nex(thisfile, 'header', hdr.orig(thischan), 'channel', 1, 'begsample', begsample, 'endsample', endsample); % always read the first and only channel
        dat(i,:) = nex.dat;
      case 'plexon_plx'
        ft_error('plx files are not supported in plexon dataset directory');
      case 'plexon_ddt'
        ft_error('ddt files are not supported in plexon dataset directory');
      otherwise
        ft_error('unsupported file in plexon dataset directory');
    end
  end
end
