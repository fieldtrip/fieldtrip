function [hdr, dat] = read_neurosim_evolution(filename, varargin)

% READ_NEUROSIM_EVOLUTION reads the "evolution" file that is written
% by Jan van der Eerden's NeuroSim software. When a directory is used
% as input, the default filename 'evolution' is read.
%
% Use as
%   [hdr, dat] = read_neurosim_evolution(filename, ...)
% where additional options should come in key-value pairs and can include
%   Vonly       = 0 or 1, only give the membrane potentials as output
%   headerOnly  = 0 or 1, only read the header information (skip the data), automatically set to 1 if nargout==1
% 
% See also FT_READ_HEADER, FT_READ_DATA

% Copyright (C) 2012 Robert Oostenveld
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
  filename = fullfile(filename, 'evolution');
end

Vonly           = ft_getopt(varargin, 'Vonly',0);
headerOnly      = ft_getopt(varargin, 'headerOnly',0);

if nargout<2 % make sure that when only one output is requested the header is returned
    headerOnly=true;
end

label = {};
orig  = {};

fid = fopen(filename, 'rb');

% read the header
line =  '#';
ishdr=1;
while ishdr==1
    % find temporal information
    if strfind(lower(line),'start time')
        dum= regexp(line, 'time\s+(\d+.\d+E[+-]\d+)', 'tokens');
        hdr.FirstTimeStamp = str2double(dum{1}{1});
    end
    if strfind(lower(line),'time bin')
        dum= regexp(line, 'bin\s+(\d+.\d+E[+-]\d+)', 'tokens');
        dt=str2double(dum{1}{1});
        hdr.Fs= 1e3/dt;
        hdr.TimeStampPerSample=dt;
    end
    if strfind(lower(line),'end time')
        dum= regexp(line, 'time\s+(\d+.\d+E[+-]\d+)', 'tokens');
        hdr.LastTimeStamp = str2double(dum{1}{1});
        hdr.nSamples=int64((hdr.LastTimeStamp-hdr.FirstTimeStamp)/dt+1);

    end
    % parse the content of the line, determine the label for each column
    colid = sscanf(line, '# column %d:', 1);
    if ~isempty(colid)
        label{colid} = [num2str(colid) rmspace(line(find(line==':'):end))];
    end
    
    offset = ftell(fid); % remember the file pointer position
    line   = fgetl(fid); % get the next line
    if ~isempty(line) && line(1)~='#' && ~isempty(str2num(line))
        % the data starts here, rewind the last line
        fseek(fid, offset, 'bof');
        line = [];
        ishdr=0;
    else
        orig{end+1} = line;
    end
end

timelab=find(~cellfun('isempty',regexp(lower(label), 'time', 'match')));

if ~headerOnly
  % read the complete data
  dat = fscanf(fid, '%f', [length(label), inf]);    
  hdr.nSamples    = length(dat(timelab, :)); %overwrites the value written in the header with the actual number of samples found
  hdr.LastTimeStamp = dat(1,end);
end

fclose(fid);

% only extract V_membrane if wanted
if Vonly
    matchLab=regexp(label,'V of (\S+) neuron','start');
    idx=find(~cellfun(@isempty,matchLab));
    if isempty(idx) % most likely a multi compartment simulation
        matchLab=regexp(label,'V\S+ of (\S+) neuron','start');
        idx=find(~cellfun(@isempty,matchLab));
    end
    if ~headerOnly
        dat=dat([timelab idx],:);
    end
    label=label([timelab idx]);
    for n=2:length(label) % renumbering of the labels
        label{n}=[num2str(n) label{n}(regexp(label{n},': V'):end)];
    end
end

% convert the header into FieldTrip style
hdr.label       = label(:);
hdr.nChans      = length(label);
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;
% also store the original ascii header details
hdr.orig        = orig(:);
[hdr.chanunit hdr.chantype] = deal(cell(length(label),1));
hdr.chantype(:) = {'evolution (neurosim)'};
hdr.chanunit(:) = {'unknown'};

function y=rmspace(x)
% remove double spaces from string
% (c) Bart Gips 2012
y=strtrim(x);
[sbeg send]=regexp(y,' \s+');
for n=1:length(sbeg)
    y(sbeg(n):send(n)-1)=[];
end

