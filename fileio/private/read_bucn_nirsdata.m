function [dat] = read_bucn_nirsdata(filename, hdr, begsample, endsample, chanindx)

% READ_BUCN_NIRSDATA reads ASCII-formatted NIRS data acquired with the
% UCL-BIRKBECK machine and postprocessed by the Paris group. The first line
% contains the channel labels and the rest of the file contains per line a
% time sample. The first column specifies the time axis.
%
% Use as
%   [dat] = read_bucn_nirsdata(filename, hdr, begsample, endsample, chanindx)
%
% See also READ_BUCN_NIRSHDR, READ_BUCN_NIRSEVENT

% Copyright (C) 2011, Jan-Mathijs Schoffelen
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

% initialize some variables
nchan = numel(hdr.label);
str   = repmat('%f', [1 nchan]);

% read the designated samples
fid = fopen(filename, 'r');
dat = textscan(fid, str, endsample-begsample+1, 'HeaderLines', begsample);
fclose(fid);

% reshape into a channelxsamples matrix
dat = cat(2, dat{:})';

% make the subselection
dat   = dat(chanindx, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the slower alternative is below which reads the whole file each time
% and thus is slow for multiple trials
%
% fid = fopen(filename, 'r');
% dat = textscan(fid, '%f', 'Headerlines', 1);
% dat = reshape(dat, nchan, []);
% dat = dat(chanindx, begsample:endsample);
%
