function [spike] = read_neurosim_spikes(filename,headerOnly)

% READ_NEUROSIM_SPIKES reads the "spikes" file that is written by Jan
% van der Eerden's NeuroSim software.  The output is represented in a
% structure that is consistent with the FieldTrip spike representation.
%
% OUTPUT
% spike: A fieldtrip raw spike structure (including header information in
% spike.hdr
% 
% INPUT
% filename: name of spike files or directory (this will default to using
% the 'spikes' file in the directory, the default neurosim naming
% convention)
% 
% headerOnly: (OPTIONAL) if this is true, only the header information is
% given directly as output, the spike data itself is not read in. (used by
% FT_READ_HEADER)
% 
% See also FT_READ_SPIKE, FT_DATATYPE_SPIKE

% Copyright (C) 2012 Robert Oostenveld, Bart Gips
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

if isdir(filename)
    filename = fullfile(filename, 'spikes');
end

if nargin<2
  headerOnly = false;
end

fid = fopen(filename, 'rb');
label = {};
orig  = {};

% read the header
line =  '#';
while ~isempty(line)
    % find temporal information
    if strfind(lower(line),'start time')
        dum= regexp(line, 'time\s+(\d+.\d+E[+-]\d+)', 'tokens');
        spike.hdr.FirstTimeStamp = str2double(dum{1}{1});
    end
    if strfind(lower(line),'time bin')
        dum= regexp(line, 'bin\s+(\d+.\d+E[+-]\d+)', 'tokens');
        dt=str2double(dum{1}{1});
        spike.hdr.Fs= 1e3/dt;
        spike.hdr.TimeStampPerSample=dt;
    end
    if strfind(lower(line),'end time')
        dum= regexp(line, 'time\s+(\d+.\d+E[+-]\d+)', 'tokens');
        spike.hdr.LastTimeStamp = str2double(dum{1}{1});
        spike.hdr.nSamples=(spike.hdr.LastTimeStamp-spike.hdr.FirstTimeStamp)/dt+1;
    end
    
    % parse the neuron labels
    colid = sscanf(line, '# %d: ', 1);
    if ~isempty(colid)
        label{colid} = strtrim(line(find(line=='#')+1:strfind(line,'position')-1));
    end
    
    offset = ftell(fid); % remember the file pointer position
    line   = fgetl(fid); % get the next line
    if ~isempty(line) && line(1)~='#'
        % the data starts here, rewind the last line
        line = [];
        fseek(fid, offset, 'bof');
    else
        orig{end+1} = line;
    end
end
% saving all lines of the original header (nothing parsed)
spike.hdr.orig=orig;

spike.hdr.nChans=length(label);
spike.hdr.nSamplesPre        = 0;
spike.hdr.nTrials            = 1;
[spike.hdr.chantype spike.hdr.chanunit] = deal(cell(length(label),1));
spike.hdr.chantype(:) = {'spike (neurosim)'};
spike.hdr.chanunit(:) = {'unknown'};

if ~headerOnly %if only the header is requested, reading in the data is not needed
% done with the header; read in all data, 
% (skipping lines starting with #, therefore, skipping the header)
dat=textscan(fid, '%f %d %s %*s %*s %*s %*s %*s', 'CommentStyle', '#');
end
fclose(fid);

if headerOnly
    % if only the header is requested the label information is written to
    % the header and the function output is just the header.
    spike.hdr.label=label;
    spike=spike.hdr;
else

% data is loaded; now put it in a structure that FieldTrip can use.
% every neuron is considered to have its own channel (identified by its
% neuron number).

spike.label=cell(1,numel(label));
spike.timestamp=cell(1,numel(label));


for n=1:numel(label)
    % write labels ([number]: [neuron type] in [network])
    spike.label{n}=strrep(label{n},'neuron ','');
    % select spike times belonging to the neuron
    sel=dat{2}==n;
    spike.timestamp{n}=dat{1}(sel)';
end
end
