function [dat] = read_bucn_nirsdata(filename, hdr, begsample, endsample, chanindx)

% READ_BUCN_NIRSDATA reads ASCII-formatted NIRS data acquired with the
% UCL-BIRKBECK machine and postprocessed by the Paris group. The first line
% contains the channel labels and the rest of the file contains per line a
% time sample. The first column specifies the time axis.
%
% Use as
%   [dat] = read_bucn_nirsdata(filename, hdr, begsample, endsample, chanindx)

% Copyright (C) 2011, Jan-Mathijs Schoffelen
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
% $Id:$

% read the whole data file apart from the header line
fid = fopen(filename, 'r');
dat = textscan(fid, '%f', 'HeaderLines', 1);
fclose(fid);

% reshape into a channelxsamples matrix
nchan = numel(hdr.label);
dat   = reshape(dat{1}, nchan, []);

% make the subselection
dat   = dat(chanindx, begsample:endsample);
