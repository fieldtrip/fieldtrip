function varargout = openvibe_mat(filename, hdr, begsample, endsample, chanindx)

% OPENVIBE_MAT reads EEG data from MATLAB file with OpenVibe data that was converted
% according to http://openvibe.inria.fr/converting-ov-files-to-matlab/
%
% Use as
%   hdr = openvibe_mat(filename);
%   dat = openvibe_mat(filename, hdr, begsample, endsample, chanindx);
%   evt = openvibe_mat(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT
% See also BIDS_TSV, BIOPAC_ACQ, BUCN_TXT, EEGSYNTH_TSV, EVENTS_TSV, LIBERTY_CSV, MAUS_TEXTGRID, MOTION_C3D, OPENBCI_TXT, OPENPOSE_KEYPOINTS, OPENSIGNALS_TXT, OPENVIBE_MAT, OPM_FIL, QUALISYS_TSV, SCCN_XDF, SENSYS_CSV, SNIRF, SPIKEGLX_BIN, UNICORN_CSV, XSENS_MVNX

% Copyright (C) 2023-2024, Robert Oostenveld
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

persistent ov previous_fullname

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

if isempty(previous_fullname) || ~isequal(fullname, previous_fullname) || isempty(ov)
  % read the header, data and events
  ov = load(fullname);
  % remember the full filename including path
  previous_fullname = fullname;
else
  % use the persistent variable to speed up subsequent read operations
end

% the MATLAB file contains 'stims', 'sampleTime', 'samples', 'samplingFreq', 'channelNames'
[nsamples, nchans] = size(ov.samples);

if needhdr
  % construct the FieldTrip header
  hdr = [];
  hdr.label       = ov.channelNames(:);
  hdr.Fs          = ov.samplingFreq;
  hdr.nChans      = nchans;
  hdr.nSamples    = nsamples;
  hdr.nTrials     = 1; % assume continuous data
  hdr.nSamplesPre = 0;
  % return the header
  varargout = {hdr};

elseif needdat
  % select the channels and samples, and transpose the result
  dat = ov.samples(begsample:endsample, chanindx)';
  % return the data
  varargout = {dat};

elseif needevt
  evt = struct();
  for i=1:size(ov.stims,1)
    evt(i).type     = 'trigger';
    evt(i).sample   = round(ov.stims(i,1)*hdr.Fs+1);
    evt(i).value    = ov.stims(i,2);
    evt(i).offset   = 0;
    evt(i).duration = 0;
  end
  % return the events
  varargout = {evt};

end % if needhdr, needdat or needevt
