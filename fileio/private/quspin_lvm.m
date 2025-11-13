function varargout = quspin_lvm(filename, hdr, begsample, endsample, chanindx)

% QUSPIN_LVM reads QuSpin OPM data from the LabVIEW measurement file.
% See also https://www.ni.com/docs/en-US/bundle/labview/page/labview-measurement-files.html
%
% Use as
%   hdr = quspin_lvm(filename);
%   dat = quspin_lvm(filename, hdr, begsample, endsample, chanindx);
%   evt = quspin_lvm(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT
% See also BIDS_TSV, BIOPAC_ACQ, BUCN_TXT, EEGSYNTH_TSV, EVENTS_TSV, LIBERTY_CSV, MAUS_TEXTGRID, MOTION_C3D, OPENBCI_TXT, OPENPOSE_KEYPOINTS, OPENSIGNALS_TXT, OPENVIBE_MAT, OPM_FIL, QUALISYS_TSV, QUSPIN_LVM, SCCN_XDF, SENSYS_CSV, SNIRF, SPIKEGLX_BIN, UNICORN_CSV, XSENS_MVNX

% Copyright (C) 2025, Robert Oostenveld
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

% determine how many header lines there are
fid = fopen(fullname, 'rt');
line = fgetl(fid);
separator = line(end); % the first line ends with the separator
headerlines = 0;
while ~startsWith(line, 'X_Value')
  headerlines = headerlines + 1;
  line = fgetl(fid);
end
fclose(fid);

if needhdr
  %% parse the header

  % read and parse the header lines
  fid = fopen(fullname, 'rt');
  for i=1:headerlines
    line = fgetl(fid);
    if startsWith(line, '***End_of_Header***')
      continue
    end
    tok = split(line, separator);
    if length(tok)==2 && ~isempty(tok{1}) && ~isempty(tok{2})
      orig.(tok{1}) = tok{2};
    elseif length(tok)>2
      sel = find(~cellfun(@isempty, tok(2:end)));
      if isscalar(sel)
        % make it a single value
        orig.(tok{1}) = tok{sel+1};
      else
        % keep it as cell-array
        orig.(tok{1}) = tok(sel+1);
      end
    end
  end

  % the first subsequent line of the file starts with X_Value and contains the channel labels
  line = fgetl(fid);
  hdr.label = split(line, separator);
  hdr.label = hdr.label(2:end-1); % drop the X_Value and Comment

  % estimate the sampling rate, we would expect this to be 375Hz or 750Hz
  sample1 = split(fgetl(fid), separator);
  sample2 = split(fgetl(fid), separator);
  t1 = str2double(sample1{1});
  t2 = str2double(sample2{1});
  hdr.Fs = 1/(t2-t1);

  % read until the end
  hdr.nSamples = 2;
  while ~feof(fid)
    fgetl(fid);
    hdr.nSamples = hdr.nSamples + 1;
  end

  % close the file
  fclose(fid);

  % construct the rest of the FieldTrip style header
  hdr.nChans = str2double(orig.Channels);
  hdr.nTrials = 1;
  hdr.nSamplesPre = 0;
  hdr.chanunit = orig.Y_Unit_Label;
  hdr.chantype = repmat({'unknown'}, 1, hdr.nChans);
  sel = startsWith(hdr.label, 'X');
  hdr.chantype(sel) = {'megmag'};
  sel = startsWith(hdr.label, 'Y');
  hdr.chantype(sel) = {'megmag'};
  sel = startsWith(hdr.label, 'Z');
  hdr.chantype(sel) = {'megmag'};
  sel = startsWith(hdr.label, 'T');
  hdr.chantype(sel) = {'digital trigger'};
  sel = startsWith(hdr.label, 'A');
  hdr.chantype(sel) = {'analog trigger'};

  % retain the original header details
  hdr.orig = orig;

  % return the header
  varargout = {hdr};

elseif needdat
  %% parse the data

  % jump to where the data starts
  fid = fopen(fullname, 'rt');
  for i=1:headerlines
    fgetl(fid);
  end
  fgetl(fid); % this is the line that starts with X_Value

  dat = zeros(length(chanindx), endsample-begsample+1);

  for i=1:(begsample-1)
    % skip these lines
    fgetl(fid);
  end

  % read the desired samples, select the desired channels
  for i=1:(endsample-begsample+1)
    sample = str2num(fgetl(fid)); %#ok<ST2NM>
    sample = sample(2:end); % the first is the timestamp
    dat(:,i) = sample(chanindx);
  end

  % close the file
  fclose(fid);

  % return the data
  varargout = {dat};

elseif needevt
  %% parse the events and trigger codes

  % use a generic function to find the flanks in the digital trigger channels
  chanindx = find(strcmp(hdr.chantype, 'digital trigger'));
  event = read_trigger(fullname, 'header', hdr, 'chanindx', chanindx);

  % return the events
  varargout = {event};

end
