function varargout = sensys_csv(filename, hdr, begsample, endsample, chanindx)

% SENSYS_CSV reads fluxgate magnetometer from the Sensys FGM3D TD system
%
% See https://sensysmagnetometer.com/products/fgm3d/
%
% Use as
%   hdr = sensys_csv(filename);
%   dat = sensys_csv(filename, hdr, begsample, endsample, chanindx);
%   evt = sensys_csv(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, MOTION_C3D, QUALISYS_TSV, LIBERTY_CSV

% Copyright (C) 2022, Robert Oostenveld
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
  fid = fopen_or_error(fullname);
  num = 1;
  header = struct();
  line = fgetl(fid);
  while ~isempty(line)
    line = fgetl(fid);
    if ~isempty(line)
      tok = strsplit(line, ',');
      header.(fixname(tok{1})) = tok{2};
    end
    num = num + 1;
  end
  fclose(fid);
  csv.header = header;
  csv.dat = readtable(fullname, 'FileType', 'text', 'NumHeaderLines', num, 'Delimiter', ',', 'ExpectedNumVariables', 5, 'VariableNamingRule', 'preserve');
else
  % use the persistent variable to speed up subsequent read operations
end

if needhdr
  %% parse the header
  hdr = [];
  hdr.Fs          = str2double(csv.header.samplerate__hz);
  hdr.nChans      = size(csv.dat,2);
  hdr.nSamples    = size(csv.dat,1);
  hdr.nSamplesPre = 1;
  hdr.nTrials     = 1;
  hdr.label       = {'timestamp', 'x', 'y', 'z', 'abs'};
  hdr.chantype    = {'misc', 'megmag', 'megmag', 'megmag', 'misc'};
  hdr.chanunit    = {'ms', 'T', 'T', 'T', 'T'};

  % return the header details
  varargout = {hdr};

elseif needdat
  [nsample, nchan] = size(csv.dat);
  dat = nan(nchan, nsample);

  %% parse the data
  % the first column contains integers, the subsequent ones data in scientific format
  % the scientific format might be either with a '.' or a ',' as decimal separator
  for i=1:nchan
    if isnumeric(csv.dat{1,i})
      dat(i,:) = table2array(csv.dat(:,i));
    else
      tmp = table2cell(csv.dat(:,i));
      tmp = strrep(tmp, ',', '.');
      dat(i,:) = str2double(tmp);
    end
  end
  % make the selection of channels and samples
  dat = dat(chanindx, begsample:endsample);

  % return the data
  varargout = {dat};

elseif needevt
  %% there are no events
  event = [];

  % return the events
  varargout = {event};

end
