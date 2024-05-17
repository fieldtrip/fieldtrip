function [varargout] = opm_fil(filename, hdr, begsample, endsample, chanindx)

% OPM_FIL reads header, data and events from OPM MEG recordings that are done at the FIL (UCL, London).
%
% Use as
%   hdr = opm_fil(filename);
%   dat = opm_fil(filename, hdr, begsample, endsample, chanindx);
%   evt = opm_fil(filename, hdr);
%
% See https://github.com/tierneytim/OPM for technical details.
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT
% See also BIDS_TSV, BIOPAC_ACQ, BUCN_TXT, EEGSYNTH_TSV, EVENTS_TSV, LIBERTY_CSV, MAUS_TEXTGRID, MOTION_C3D, OPENBCI_TXT, OPENPOSE_KEYPOINTS, OPENSIGNALS_TXT, OPENVIBE_MAT, OPM_FIL, QUALISYS_TSV, SCCN_XDF, SENSYS_CSV, SNIRF, SPIKEGLX_BIN, UNICORN_CSV, XSENS_MVNX

% Copyright (C) 2020-2024, Robert Oostenveld
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
needevt = (nargin==2);
needdat = (nargin==5);

% It is permitted to pass any of the files that comprise the dataset
% figure out which files comprise the dataset
[p, f, x] = fileparts(filename);
assert(isfile(filename), sprintf('file "%s" not found', filename));

% Remove the _meg part of the filename
if endsWith(f, '_meg')
  f = f(1:end-4);
end

% Get the BIDS compliant files
channelsfile  = fullfile(p,[f '_channels.tsv']);
coordsysfile  = fullfile(p,[f '_coordsystem.json']);
datafile      = fullfile(p,[f '_meg.bin']);
headerfile    = fullfile(p,[f '_meg.json']);
positionsfile = fullfile(p,[f '_positions.tsv']);

if needhdr
  %% read the header
  try
    fid = fopen(headerfile, 'rt');
    header = jsondecode(fread(fid, [1 inf], 'char=>char'));
    fclose(fid);
  catch
    ft_warning('Cannot open header');
  end

  % determine the precision from the json header
  if isfield(header, 'Precision')
    precision = header.Precision;
  else
    precision = 'single';
  end

  switch precision
    case 'single'
      samplesize = 4;
    case 'double'
      samplesize = 8;
  end

  channels  = readtable(channelsfile, 'Delimiter', 'tab', 'FileType', 'text');

  % this one is optional
  if exist(coordsysfile, 'file')
    fid = fopen(coordsysfile, 'rt');
    coordsys = jsondecode(fread(fid, [1 inf], 'char=>char'));
    fclose(fid);
  else
    coordsys = [];
  end

  % this one is optional
  if exist(positionsfile, 'file')
    positions = readtable(positionsfile, 'Delimiter', 'tab', 'FileType', 'text');
  else
    positions = [];
  end

  d = dir(datafile);

  hdr.label       = channels.name;
  hdr.nChans      = size(channels, 1);
  hdr.nSamples    = (d.bytes/(size(channels, 1) * samplesize))-1;
  hdr.nSamplesPre = 0; % continuous data
  hdr.nTrials     = 1; % continuous data
  hdr.Fs          = header.SamplingFrequency;
  hdr.chantype    = channels.type;
  try
    hdr.chanunit    = channels.unit;
  catch
    hdr.chanunit    = repmat({'unknown'}, size(hdr.label));
  end

  % keep the original header details
  hdr.orig.header     = header;
  hdr.orig.coordsys   = coordsys;
  hdr.orig.positions  = positions;
  hdr.orig.channels   = channels;

  if ~isempty(positions)
    % FIXME construct a grad structure
    hdr.grad.label = positions.name;
    hdr.grad.coilpos = [positions.Px positions.Py positions.Pz];
    hdr.grad.coilori = [positions.Ox positions.Oy positions.Oz];
    hdr.grad.type = 'meg';
    hdr.grad.tra  = diag(ones(1,length(positions.name)));
    if ~isempty(coordsys)
      hdr.grad.unit = coordsys.MEGCoordinateUnits;
      hdr.grad.coordsys = coordsys.MEGCoordinateSystem;
    end
  end

  % return the header details
  varargout = {hdr};

elseif needdat
  %% read the data, note that it is big endian

  % determine the precision from the json header
  % see https://github.com/fieldtrip/fieldtrip/issues/2240
  if isfield(hdr.orig.header, 'Precision')
    precision = hdr.orig.header.Precision;
  else
    precision = 'single';
  end

  switch precision
    case 'single'
      samplesize = 4;
    case 'double'
      samplesize = 8;
  end

  fid = fopen(datafile, 'rb');
  fseek(fid, begsample*samplesize*hdr.nChans, 'bof');
  dat = fread(fid,[hdr.nChans, (endsample-begsample+1)], precision, 0, 'b');
  fclose(fid);

  % take the selection of channels
  dat = dat(chanindx,:);

  % return the data
  varargout = {dat};

elseif needevt
  %% read the events

  % FIXME this does not yet allow the user to override the defaults for ft_read_event, such as detectflank
  event = read_trigger(filename, 'header', hdr, 'dataformat', 'opm_fil', 'chanindx', find(strcmpi(hdr.chantype, 'TRIG')));

  % return the events
  varargout = {event};

end
