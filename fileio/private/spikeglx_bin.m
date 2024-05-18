function varargout = spikeglx_bin(filename, hdr, begsample, endsample, chanindx)

% SPIKEGLX_BIN reads Neuropixel data from SpikeGLX .bin files
%
% See https://github.com/jenniferColonell/SpikeGLX_Datafile_Tools
%
% Use as
%   hdr = spikeglx_bin(filename);
%   dat = spikeglx_bin(filename, hdr, begsample, endsample, chanindx);
%   evt = spikeglx_bin(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT
% See also BIDS_TSV, BIOPAC_ACQ, BUCN_TXT, EEGSYNTH_TSV, EVENTS_TSV, LIBERTY_CSV, MAUS_TEXTGRID, MOTION_C3D, OPENBCI_TXT, OPENPOSE_KEYPOINTS, OPENSIGNALS_TXT, OPENVIBE_MAT, OPM_FIL, QUALISYS_TSV, SCCN_XDF, SENSYS_CSV, SNIRF, SPIKEGLX_BIN, UNICORN_CSV, XSENS_MVNX

% Copyright (C) 2024, Robert Oostenveld
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

% ensure that the external toolbox is available
ft_hastoolbox('spikeglx', 1);

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

% use the full filename including path to distinguish between similarly named files in different directories
[p, f, x] = fileparts(filename);

if strcmp(x, '.meta')
  x = '.bin';
end

if isempty(p)
  % no path was specified
  fullname = which(filename);
elseif startsWith(p, ['.' filesep])
  % a relative path was specified
  fullname = fullfile(pwd, p(3:end), [f, x]);
else
  fullname = filename;
end


if needhdr
  [p, f, x] = fileparts(fullname);
  meta = SGLX_readMeta.ReadMeta([f x], p);
  chan = split(meta.snsChanMap, ')(');

  hdr = [];
  hdr.nChans = str2double(meta.nSavedChans);
  hdr.Fs = SGLX_readMeta.SampRate(meta);
  hdr.nSamples = str2double(meta.fileSizeBytes)/(hdr.nChans*2);
  hdr.nSamplesPre = 0;  % continuous data
  hdr.nTrials = 1;      % continuous data
  hdr.label = strtok(chan(2:end), ';'); % the first does not seem to be a channel
  hdr.orig = meta;      % keep the original details

  % return the header
  varargout = {hdr};

elseif needdat
  % the header has already been read and passed as input argument
  meta = hdr.orig;

  % read the binary data
  dat = SGLX_readMeta.ReadBin(begsample-1, endsample-begsample+1, meta, [f x], p);

  % apply the calibration
  if strcmp(meta.typeThis, 'imec')
    dat = SGLX_readMeta.GainCorrectIM(dat, 1:size(dat,1), meta);
  else
    dat = SGLX_readMeta.GainCorrectNI(dat, 1:size(dat,1), meta);
  end

  % make the selection of channels
  dat = dat(chanindx,:);

  % return the data
  varargout = {dat};

elseif needevt
  % the header has already been read and passed as input argument
  meta = hdr.orig;

  % For a digital channel: read this digital word dw in the saved file
  % (1-based). For imec data there is never more than one saved digital word.
  dw = 1;

  ntrig = 8;

  dat = SGLX_readMeta.ReadBin(0, hdr.nSamples, meta, [f x], p);
  dig = SGLX_readMeta.ExtractDigital(dat, meta, dw, 1:ntrig);

  event = [];
  for i=1:ntrig
    trig = single(dig(i,:));
    onset  = find(diff([0 trig])>0);
    offset = find(diff([trig 0])<0);
    for j=1:numel(onset)
      event(end+1).type = 'trigger';
      event(end  ).value = i;
      event(end  ).sample = onset(j);
      event(end  ).duration = offset(j)-onset(j)+1;
      event(end  ).offset = 0;
    end
  end

  % return the events
  varargout = {event};
end