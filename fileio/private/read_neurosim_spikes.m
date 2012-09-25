function [spike] = read_neurosim_spikes(filename)

% READ_NEUROSIM_SPIKES reads the "spikes" file that is written by Jan
% van der Eerden's NeuroSim software.  The output is represented in a
% structure that is consistent with the FieldTrip spike representation.
%
% See also FT_READ_SPIKE, FT_DATATYPE_SPIKE

% Copyright (C) 2012 Robert Oostenveld
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

% count the number of lines in the file
numlines = 0;
num = 1024;
while num==1024
  [str, num] = fread(fid, num, 'uchar');
  numlines = numlines + sum(str==10);
end
fseek(fid, 0, 'bof');

line = '#';
while (line(1)=='#')
  prevline = line;
  line = fgetl(fid);
  numlines = numlines - 1; % reduce the number of data lines
end
numlines = numlines + 1; % the line that was just read actually contains data
numlines = numlines - 1; % the last line of the file is usually an empty line, or is not properly defined

% try to parse the 
header  = tokenize(prevline, ' ', 1);
header  = header(2:end); % the first element is '#'

% determine the number of data columns
tok     = tokenize(line, ' ', 1);
islabel = cellfun(@isempty, cellfun(@str2num, tok, 'UniformOutput', false));

numeric = zeros(numlines,sum(~islabel));
label   = cell(numlines,sum(islabel));

ft_progress('init', 'etf', 'reading ascii file...');
for i=1:numlines
  ft_progress(i/numlines, 'reading line %d from %d', i, numlines)
  tok = tokenize(line, ' ', 1);
  numeric(i,:) = str2double(tok(~islabel));
  label(i,:)   = tok(islabel);
  line = fgetl(fid);
end
ft_progress('close');
fclose(fid);

% Having all the data, the next is to get it all represented according
% to the FieldTrip spike representation, see T_DATATYPE_SPIKE

% A minimal representation would be to pretend that all neurons were measured
% with a single electrode/channel and that they were afterwards sorted in the
% different units
%
% spike.label     = {'spikes'}
% spike.timestamp = {numeric(:,1)'};
% spike.unit      = {numeric(:,2)'};

% it is more convenient to represent each neuron in its own channel
number          = unique(numeric(:,2)');
spike.label     = cell(1, length(number), 1);
spike.timestamp = cell(1, length(number), 1);
for i=1:length(number)
  sel = numeric(:,2)==number(i);
  spike.label{i}     = sprintf('%d', number(i));
  spike.timestamp{i} = numeric(sel,1)';
end

