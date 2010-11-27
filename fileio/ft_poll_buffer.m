function numbers = ft_poll_buffer(filename, thresholds, timeout)
% FT_POLL_BUFFER polls the FieldTrip buffer described by 'filename' for
% the number of events and samples. You can also use this function to
% wait up to a given amount of milliseconds until the number of samples
% OR the number of events present in the buffer is bigger than a specified
% threshold. The output of this function is a struct with elements
% 'nsamples' and 'nevents', that always contain the current quantities.
% (So in case of timeouts, these might be smaller than the thresholds).
%
% Examples
% --------
% The following will immediately return the two quantities:
% >> numbers = ft_poll_buffer('buffer://localhost:1972');
%
% The following will wait up to 50 milliseconds for more than 100 samples in the buffer.
%   thr = []; thr.nsamples = 100; timeout = 50;
%   numbers = ft_poll_buffer('buffer://localhost:1972', thr, timeout);
%
% The following will wait up to 80 milliseconds for more than 30 events in the buffer.
%   thr = []; thr.nevents = 30; timeout = 80;
%   numbers = ft_poll_buffer('buffer://localhost:1972', thr, timeout);
%
% Repeatedly wait for new data (either samples or events):
%
%   lastNum = []; 
%   lastNum.nevents = 0; 
%   lastNum.nsamples = 0; 
%   timeout = 100;
%   while someCriterion
%      newNum = ft_poll_buffer('buffer://localhost:1972', lastNum, timeout);
%      if newNum.nsamples > lastNum.nsamples
%         ... read and process samples ...
%      if newNum.nevents > lastNum.nevents
%         ... read and process events ...
%      lastNum = newNum;
%   end

% Copyright (C) 2010, Stefan Klanke
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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