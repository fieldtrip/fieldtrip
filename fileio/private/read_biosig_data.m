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
% $Log: read_biosig_data.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.1  2004/03/10 11:50:31  roberto
% new implementation of wrapper functions around the BIOSIG toolbox (specifically section T200)
%

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
