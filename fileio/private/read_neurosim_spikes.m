function [spike] = read_neurosim_spikes(filename)

% READ_NEUROSIM_SPIKES reads the "spikes" file that is written by Jan
% van der Eerden's NeuroSim software.  The output is represented in a
% structure that is consistent with the FieldTrip spike representation.
%
% See also FT_READ_SPIKE, FT_DATATYPE_SPIKE

% Copyright (C) 2012 Robert Oostenveld, Bart Gips
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

if isdir(filename)
    filename = fullfile(filename, 'spikes');
end
fid = fopen(filename, 'rb');


dat=textscan(fid, '%f %d %s %*s %*s %*s %*s %*s', 'CommentStyle', '#');
fclose(fid);

% data is loaded; now put it in a structure that FieldTrip can use.
% every neuron is considered to have its own channel (identified by its
% neuron number).

[number, ~, idx]=unique(dat{2});
spike.label=cell(1,length(number));
spike.timestamp=cell(1,length(number));


for n=1:length(number)
    spike.label{n}=num2str(number(n));
    sel=idx==number(n);
    spike.timestamp{n}=dat{1}(sel)';
end
