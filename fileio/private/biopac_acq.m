function varargout = biopac_acq(filename, hdr, begsample, endsample, chanindx)

% BIOPAC_ACQ is a wrapper to for the reading function from Mathworks file exchange.
%
% Use as
%   hdr = biopac_acq(filename);
%   dat = biopac_acq(filename, hdr, begsample, endsample, chanindx);
%   evt = biopac_acq(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT

% Copyright (C) 2018 Robert Oostenveld
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

persistent acq previous_fullname

% add fieldtrip/external/fileexchange to the path
ft_hastoolbox('fileexchange', 1);

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

% use the full filename including path to distinguish between similarly named files in different directories
fullname = which(filename);

if isempty(previous_fullname) || ~isequal(fullname, previous_fullname)
  % remember the full filename including path
  previous_fullname = fullname;
  % read the header and data
  acq = load_acq(filename, false);
else
  % use the persistent variable to speed up subsequent read operations
end

if needhdr
  % convert to FieldTrip header representation
  hdr.Fs           = 1000/acq.hdr.graph.sample_time;
  hdr.nChans       = size(acq.data,2);
  hdr.nSamples     = size(acq.data,1);
  hdr.nSamplesPre  = 0;
  hdr.nTrials      = 1; % assume that it is a single continuous segment
  for i=1:numel(acq.hdr.per_chan_data)
    hdr.label{i} = acq.hdr.per_chan_data(i).comment_text;
  end
  hdr.orig = acq.hdr;
  varargout{1} = hdr;
  
elseif needevt
  % convert to FieldTrip event representation
  event = [];
  varargout{1} = event;
  
elseif needdat
  % select the requested channels and samples from the data and transpose
  dat = acq.data(begsample:endsample,chanindx)';
  varargout{1} = dat;
end
