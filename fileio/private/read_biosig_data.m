function [dat] = read_biosig_data(filename, hdr, begsample, endsample, chanindx)

% READ_BIOSIG_DATA reads data from EEG file using the BIOSIG
% toolbox and returns it in the FCDC framework standard format
%
% Use as
%  [dat] = read_biosig_data(filename, hdr, begsample, endsample, chanindx)
% where the header has to be read before with READ_BIOSIG_HEADER.
%
% The following data formats are supported: EDF, BKR, CNT, BDF, GDF

% Copyright (C) 2004-2012, Robert Oostenveld
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

persistent cacheheader        % for caching

% get the details from the file
details = dir(filename);
if nargin>4
  % the channel selection is made when opening the file
  % and should be remembered over multiple calls
  details.chan = chanindx;
end

% open the file and read the header or use the cached header
if isempty(cacheheader) || ~isequal(details, cacheheader.details)
  % close previous file, if different
  if ~isempty(cacheheader) && exist(cacheheader.fullname, 'file')
    sclose(cacheheader);
  end
  if nargin>4
    HDR = sopen(filename,'r', chanindx);
  else
    HDR = sopen(filename,'r');
  end
  % put the header in the cache
  cacheheader = HDR;
  % update the header details (including time stampp, size and name)
  cacheheader.details = dir(filename);
  % we need full path for file closing
  cacheheader.fullname = filename;
else
  HDR = cacheheader;
end

begtime = (begsample-1) / hdr.Fs;
endtime = (endsample  ) / hdr.Fs;
duration = endtime-begtime;

% read the selected data
dat = sread(HDR, duration, begtime);
dat = dat';
