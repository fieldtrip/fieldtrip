function [dat] = read_neuralynx_cds(filename, hdr, begsample, endsample, chanindx);

% READ_NEURALYNX_CDS reads selected samples and channels from a  combined Neuralynx dataset with separate subdirectories for the LFP, MUA and spike channels
%
% Use as
%    hdr = read_neuralynx_cds(parentdir)
%    dat = read_neuralynx_cds(parentdir, hdr, begsample, endsample, chanindx)

% Copyright (C) 2006, Robert Oostenveld
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

needhdr = (nargin==1);
needdat = (nargin>=2);

dirlist = dir(filename);
haslfp   = any(filetype_check_extension({dirlist.name}, 'lfp'));
hasmua   = any(filetype_check_extension({dirlist.name}, 'mua'));
hasspike = any(filetype_check_extension({dirlist.name}, 'spike'));
%hasttl  = any(filetype_check_extension({dirlist.name}, 'ttl'));   % separate file with original Parallel_in
%hastsl  = any(filetype_check_extension({dirlist.name}, 'tsl'));   % separate file with original TimeStampLow
%hastsh  = any(filetype_check_extension({dirlist.name}, 'tsh'));   % separate file with original TimeStampHi

if needhdr
  % read the header from the LFP dataset
  if haslfp
    lfpdir = fullfile(filename, dirlist(find(filetype_check_extension({dirlist.name}, 'lfp'))).name);
    lfphdr = read_header(lfpdir);
    lfpfil = lfphdr.filename;
    for i=1:lfphdr.nChans
      lfplab{i} = sprintf('lfp_%s', lfphdr.label{i});
    end
  else
    lfphdr = [];
    lfpfil = {};
    lfplab = {};
  end

  % read the header from the MUA dataset
  if hasmua
    muadir = fullfile(filename, dirlist(find(filetype_check_extension({dirlist.name}, 'mua'))).name);
    muahdr = read_header(muadir);
    muafil = muahdr.filename;
    for i=1:muahdr.nChans
      mualab{i} = sprintf('mua_%s', muahdr.label{i});
    end
  else
    muahdr = [];
    muafil = {};
    mualab = {};
  end

  % read the header from the SPIKE dataset
  if hasspike
    spikedir = fullfile(filename, dirlist(find(filetype_check_extension({dirlist.name}, 'spike'))).name);
    spikehdr = read_header(spikedir);
    spikefil = spikehdr.filename;
    for i=1:spikehdr.nChans
      spikelab{i} = sprintf('spike_%s', spikehdr.label{i});
    end
  else
    spikehdr = [];
    spikefil = {};
    spikelab = {};
  end

  if haslfp
    % assume that the LFP header describes the MUA and the spikes as well
    hdr = lfphdr;
  else hasmua
    % assume that the MUA header describes the spikes as well
    hdr = muahdr;
  end

  % concatenate the lfp/mua/spike channels, they have already been renamed
  hdr.label    = cat(1, lfplab(:), mualab(:), spikelab(:));
  hdr.filename = cat(1, lfpfil(:), muafil(:), spikefil(:));
  hdr.nChans = length(hdr.label);
  % remember the original header details
  hdr.orig    = {};
  hdr.orig{1} = lfphdr;
  hdr.orig{2} = muahdr;
  hdr.orig{3} = spikehdr;

  % return the header
  dat = hdr;
else

  % use the original header details
  lfphdr   = hdr.orig{1};
  muahdr   = hdr.orig{2};
  spikehdr = hdr.orig{3};
  % the spike header does not contain these, but requires them from the LFP or MUA
  spikehdr.FirstTimeStamp     = hdr.FirstTimeStamp;
  spikehdr.TimeStampPerSample = hdr.TimeStampPerSample;

  if ~isempty(lfphdr)
    Nlfp = lfphdr.nChans;
  else
    Nlfp = 0;
  end

  if ~isempty(muahdr)
    Nmua = muahdr.nChans;
  else
    Nmua = 0;
  end

  if ~isempty(spikehdr)
    Nspike = spikehdr.nChans;
  else
    Nspike = 0;
  end

  % determine which files should be read from each respective dataset
  selchan = zeros(hdr.nChans,1);
  selchan(chanindx) = 1;
  sellfp   = find(selchan(1:Nlfp));
  selmua   = find(selchan((Nlfp+1):(Nlfp+Nmua)));
  selspike = find(selchan((Nlfp+Nmua+1):(Nlfp+Nmua+Nspike)));

  if haslfp 
    % the header contains the relevant file names
    lfpdat = read_neuralynx_ds([], lfphdr, begsample, endsample, sellfp);
  else
    lfpdat = [];
  end

  if hasmua  
    % the header contains the relevant file names
    muadat = read_neuralynx_ds([], muahdr, begsample, endsample, selmua);
  else
    muadat = [];
  end

  if hasspike  
    % the header contains the relevant file names
    spikedat = read_neuralynx_ds([], spikehdr, begsample, endsample, selspike);
  else
    spikedat = [];
  end

  % combine all data into a single array
  dat = cat(1, lfpdat, muadat, spikedat);

end

