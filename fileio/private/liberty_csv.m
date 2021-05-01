function varargout = liberty_csv(filename, hdr, begsample, endsample, chanindx)

% LIBERTY_CSV reads motion capture data from the Polhemus Liberty system
%
% Use as
%   hdr = liberty_csv(filename);
%   dat = liberty_csv(filename, hdr, begsample, endsample, chanindx);
%   evt = liberty_csv(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, MOTION_C3D, QUALISYS_TSV

% Copyright (C) 2021, Robert Oostenveld
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

persistent csv previous_fullname

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

% use the full filename including path to distinguish between similarly named files in different directories
[p, f, x] = fileparts(filename);
if isempty(p)
  % no path was specified
  fullname = which(filename);
elseif startsWith(p, ['.' filesep])
  % a relative path was specified
  fullname = fullfile(pwd, p(3:end), [f, x]);
else
  fullname = filename;
end

if isempty(previous_fullname) || ~isequal(fullname, previous_fullname) || isempty(csv)
  % remember the full filename including path
  previous_fullname = fullname;
  % read and remember the file content
  csv = readtable(fullname, 'FileType', 'text', 'Delimiter', ',', 'ReadVariableNames', true);
  csv.Properties.VariableNames = cellfun(@lower, csv.Properties.VariableNames, 'UniformOutput', false);
else
  % use the persistent variable to speed up subsequent read operations
end

nsensor = max(csv.sensor);
selpos = find(startsWith(csv.Properties.VariableNames, 'position'));
selori = find(startsWith(csv.Properties.VariableNames, 'orientation'));
selother = setdiff(1:length(csv.Properties.VariableNames), [selpos selori]);

if needhdr
  %% parse the header
  
  hdr = [];
  hdr.nChans = length(csv.Properties.VariableNames) - 6 + 6*nsensor; % three per sensor position, three per sensor orientation
  hdr.nSamples = csv.sample(end) - csv.sample(1) + 1;
  hdr.nSamplesPre = 0;    % assume continuous data
  hdr.nTrials = 1;        % assume continuous data
  
  hdr.label = csv.Properties.VariableNames(selother)';
  offset = length(selother);
  for i=1:nsensor
    hdr.label{offset + (i-1)*6+1} = sprintf('position_%d_x', i);
    hdr.label{offset + (i-1)*6+2} = sprintf('position_%d_y', i);
    hdr.label{offset + (i-1)*6+3} = sprintf('position_%d_z', i);
    hdr.label{offset + (i-1)*6+4} = sprintf('orientation_%d_x', i);
    hdr.label{offset + (i-1)*6+5} = sprintf('orientation_%d_y', i);
    hdr.label{offset + (i-1)*6+6} = sprintf('orientation_%d_z', i);
  end
  
  timestamp = mean(reshape(csv.timestamp, nsensor, []), 1);
  p = polyfit(0:length(timestamp)-1, timestamp, 1);
  % p(1) is the slope in milliseconds
  % p(2) is the offset
  hdr.Fs = 1000 / p(1);
  
  hdr.chantype = repmat({'unknown'}, size(hdr.label));
  hdr.chanunit = repmat({'unknown'}, size(hdr.label));
  hdr.chantype(ismember(hdr.label, {'marker', 'sync'})) = {'trigger'};
  hdr.chantype(startsWith(hdr.label, {'position'})) = {'position'};
  hdr.chantype(startsWith(hdr.label, {'orientation'})) = {'orientation'};
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% parse the data
  dat = nan(hdr.nChans, hdr.nSamples);
  for i=1:length(selother)
    tmp = csv{:,i};
    tmp = reshape(tmp, nsensor, []);
    tmp = mean(tmp, 1);
    dat(i,:) = tmp;
  end
  
  offset = length(selother);
  
  % reshape the positions and insert them in the data matrix
  tmp = csv{:,selpos};
  tmp = reshape(tmp', 3*nsensor, []);
  chansel = (6*repmat(1:nsensor, 3, 1)-6) + repmat(1:3, nsensor, 1)' + offset;
  dat(chansel,:) = tmp;
  
  % reshape the orientations and insert them in the data matrix
  tmp = csv{:,selori};
  tmp = reshape(tmp', 3*nsensor, []);
  chansel = (6*repmat(1:nsensor, 3, 1)-3) + repmat(1:3, nsensor, 1)' + offset;
  dat(chansel,:) = tmp;
  
  % make a selection of channels and samples
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% parse the events, this will use a recursive call to get the trigger channel data
  trigindx = find(ismember(hdr.label, {'marker', 'sync'}));
  event = read_trigger(fullname, 'chanindx', trigindx, 'detectflank', 'both');
  
  % return the events
  varargout = {event};
  
end
