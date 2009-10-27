function [dat] = read_biosig_data(filename, hdr, begsample, endsample, chanindx);

% READ_BIOSIG_DATA reads data from EEG file using the BIOSIG
% toolbox and returns it in the FCDC framework standard format
%
% Use as
%  [dat] = read_biosig_data(filename, hdr, begsample, endsample, chanindx)
% where the header has to be read before with READ_BIOSIG_HEADER.
%
% The following data formats are supported: EDF, BKR, CNT, BDF, GDF

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Copyright (C) 2004, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% open the file, read the header
% it should be the same as hdr.orig
biosig = sopen(filename,'r');

begtime = (begsample-1) / hdr.Fs;
endtime = (endsample  ) / hdr.Fs;
duration = endtime-begtime;

% read the selected data and close the file
dat = sread(biosig, duration, begtime);
dat = dat';
sclose(biosig);

if nargin>4
  % select the channels of interest
  dat = dat(chanindx,:);
end
