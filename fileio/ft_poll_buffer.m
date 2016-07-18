function numbers = ft_poll_buffer(filename, thresholds, timeout)

% FT_POLL_BUFFER is deprecated.
%
% Please use FT_READ_DATA and FT_READ_EVENT with the  'blocking' and
% the 'timeout' options.

% Copyright (C) 2010, Stefan Klanke
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

threshtime = [-1 -1 0];
if nargin>1
	if isfield(thresholds, 'nsamples')
		threshtime(1) = thresholds.nsamples;
	end
	if isfield(thresholds, 'nevents')
		threshtime(2) = thresholds.nevents;
	end
end
if nargin>2
    threshtime(3) = timeout;
end

% determine the filetype
dataformat = ft_filetype(filename);
switch dataformat
  case 'fcdc_buffer'
    % read from a networked buffer for realtime analysis
    [host, port] = filetype_check_uri(filename);
    numbers = buffer('wait_dat', threshtime, host, port);
  otherwise
    error 'Polling not supported for this data format.';
end
