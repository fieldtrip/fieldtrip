function [event] = read_bucn_nirsevent(filename)

% READ_BUCN_NIRSEVENT reads the event information of ASCII-formatted NIRS 
% data acquired with the UCL-BIRKBECK machine and postprocessed by the
% Paris group. The first line contains the header-info and the rest of
% the file contains per line an event. The first column specifies the
% time of the event in samples, the second column specifies the time of the
% event in seconds, the third column contains the event type and the fourth
% column is the event value.
%
% Use as
%   [event] = read_bucn_nirshdr(filename)
%
% See also READ_BUCN_NIRSHDR, READ_BUCN_NIRSDATA

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

fid = fopen_or_error(filename, 'r');

% read the first line
line1 = textscan(fid,         '%[^\n]',1);
label = textscan(line1{1}{1}, '%[^\t]');
label = label{1};
ncol  = numel(label);
str   = ['%f%f', repmat('%s', [1 ncol-2])];

% read the rest
dat = textscan(fid, str);
fclose(fid);

% sample numbers seem to be 0-based, FieldTrip uses matlab convention, i.e.
% 1-based
dat{1} = dat{1}+1;

% convert the numeric array into a cell-array to facilitate storing it in
% a struct-array
nevent = numel(dat{1});
dat{1} = mat2cell(dat{1}, ones(1,nevent), 1);

% create struct-array
event  = struct('sample', dat{1}, 'type', dat{3}, 'value', dat{4});
