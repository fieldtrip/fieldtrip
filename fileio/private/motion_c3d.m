function varargout = motion_c3d(filename, hdr, begsample, endsample, chanindx)

% MOTION_C3D reads motion tracking data from a C3D file, see https://www.c3d.org
%
% Use as
%   hdr = motion_c3d(filename);
%   dat = motion_c3d(filename, hdr, begsample, endsample, chanindx);
%   evt = motion_c3d(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV

% Copyright (C) 2019 Robert Oostenveld
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

persistent c3d previous_fullname

% ensure that the mex files from the external C3D toolbox are available, see https://github.com/pyomeca/ezc3d
ft_hastoolbox('ezc3d', 1);

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

if isempty(previous_fullname) || ~isequal(fullname, previous_fullname) || isempty(c3d)
  % remember the full filename including path
  previous_fullname = fullname;
  % read the header and data
  c3d = ezc3dRead(fullname);
else
  % use the persistent variable to speed up subsequent read operations
end

if needhdr
  %% parse the header
  
  % In release version 1.1.0 of ezc3d the c3d.parameters.POINT.LABELS field
  % was a cell-array. In version 1.2.4 it became a structure, that contains
  % a subfield DATA as a cell-array. Also the numeric fields in 1.2.4 have the DATA
  % subfield.
  
  hdr = [];
  hdr.label = {};
  for i=1:numel(c3d.parameters.POINT.LABELS.DATA)
    hdr.label{end+1} = [c3d.parameters.POINT.LABELS.DATA{i} '_x'];
    hdr.label{end+1} = [c3d.parameters.POINT.LABELS.DATA{i} '_y'];
    hdr.label{end+1} = [c3d.parameters.POINT.LABELS.DATA{i} '_z'];
  end
  hdr.nChans      = numel(hdr.label);
  hdr.nSamples    = c3d.parameters.POINT.FRAMES.DATA;
  hdr.nSamplesPre = 0; % continuous data
  hdr.nTrials     = 1; % continuous data
  hdr.Fs          = c3d.parameters.POINT.RATE.DATA;
  hdr.chantype    = repmat({'motion'}, size(hdr.label));
  hdr.chanunit    = repmat(c3d.parameters.POINT.UNITS.DATA, size(hdr.label));
  
  % keep a copy of the original header details
  hdr.orig = c3d.parameters;
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% parse the data
  
  % this only deals with one type of data, the file format itself is much
  % more complex
  
  nchan   = numel(c3d.parameters.POINT.LABELS.DATA)*3;
  nsample = c3d.parameters.POINT.FRAMES.DATA;
  dat = reshape(c3d.data.points, [nchan nsample]);
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% parse the events
  ft_warning('reading of events is not yet implemented');
  
  % return the events
  varargout = {[]};
  
end
