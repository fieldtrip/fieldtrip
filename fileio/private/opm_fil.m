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
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D, EVENTS_TSV

% Copyright (C) 2020, Robert Oostenveld
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

channelsfile  = fullfile(p, 'channels.tsv');
coordsysfile  = fullfile(p, 'coordsystem.json');
datafile      = fullfile(p, 'meg.bin');
headerfile    = fullfile(p, 'meg.json');
positionsfile = fullfile(p, 'positions.tsv');

precision = 'double';
switch precision
  case 'single'
    samplesize = 4;
  case 'double'
    samplesize = 8;
end


if needhdr
  %% read the header
  
  fid = fopen(headerfile, 'rt');
  header = jsondecode(fread(fid, [1 inf], 'char=>char'));
  fclose(fid);
  
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
  hdr.nSamples    = d.bytes/(size(channels, 1) * samplesize);
  hdr.nSamplesPre = 0; % continuous data
  hdr.nTrials     = 1; % continuous data
  hdr.Fs          = header.SamplingFrequency;
  hdr.chantype    = channels.type;
  hdr.chanunit    = repmat({'unknown'}, size(hdr.label));
  
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
    if ~isempty(coordsys)
      hdr.grad.unit = coordsys.MEGCoordinateUnits;
      hdr.grad.coordsys = coordsys.MEGCoordinateSystem;
    end
  end
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% read the data, note that it is big endian
  
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
