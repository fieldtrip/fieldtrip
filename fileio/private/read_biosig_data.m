function [dat] = read_biosig_data(filename, hdr, begsample, endsample, chanindx);

% READ_BIOSIG_DATA reads data from EEG file using the BIOSIG
% toolbox and returns it in the FCDC framework standard format
%
% Use as
%  [dat] = read_biosig_data(filename, hdr, begsample, endsample, chanindx)
% where the header has to be read before with READ_BIOSIG_HEADER.
%
% The following data formats are supported: EDF, BKR, CNT, BDF, GDF

% Copyright (C) 2004, Robert Oostenveld
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

persistent cacheheader        % for caching

% open the file and read the header or use cached header
details = dir(filename);
if isempty(cacheheader) || ~isequal(details, cacheheader.details)
    % close previous file, if different
    if ~isempty(cacheheader) && exist(cacheheader.fullname, 'file'), sclose(cacheheader); end 
    biosig = sopen(filename,'r'); 
    % put the header in the cache
    cacheheader = biosig;
    % update the header details (including time stampp, size and name)
    cacheheader.details = dir(filename);
    % we need full path for file closing
    cacheheader.fullname = filename;
else
    biosig = cacheheader;  
end

begtime = (begsample-1) / hdr.Fs;
endtime = (endsample  ) / hdr.Fs;
duration = endtime-begtime;

% read the selected data 
dat = sread(biosig, duration, begtime);
dat = dat';

if nargin>4
  % select the channels of interest
  dat = dat(chanindx,:);
end
