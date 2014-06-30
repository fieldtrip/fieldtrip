function [swf] = read_besa_swf(filename)

% READ_BESA_SWF
%
% Use as
%   [swf] = read_besa_swf(filename)
%
% This will return a structure with the header information in
%   swf.label     cell array with labels
%   swf.data      data matrix, Nchan X Npnts
%   swf.npnt
%   swf.tsb
%   swf.di
%   swf.sb

% Copyright (C) 2006, Robert Oostenveld
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

fid = fopen(filename);
line = fgetl(fid);
head = sscanf(line,'Npts= %f TSB= %f DI= %f SB= %f'); % new BESA versions include more information in the .swf file header; skip that
Npts = head(1);
TSB  = head(2);
DI   = head(3);
SB   = head(4);

offset = ftell(fid);

% try reading a label and a number
dum = textscan(fid, '%s %f', 1);
fseek(fid, offset, 'bof');

if  isempty(dum{2})
  % first line contains all labels, each subsequent line is one timepoint for all channels
  line = fgetl(fid);
  line(line==' ') = [];
  label = tokenize(line, ':')';
  if isempty(label{end})
    label = label(1:(end-1));  % the last one is empty
  end
  Nsrc = length(label);
  dat = fscanf(fid, '%f', inf);
  dat = reshape(dat, Nsrc, Npts);
else
  % each line starts with a single label, followed by all timepoints for that channel
  dat = [];
  label = {};
  count = 0;
  while ~feof(fid)
    count = count+1;
    % read one line at a time
    line = fgetl(fid);
    sel = find(line==':', 1);
    label{count} = deblank(line(1:(sel-1)));
    buf = sscanf(line((sel+1):end), '%f');
    dat(count,:) = buf(:);
  end
end

fclose(fid);

% assign the output, should be similar as READ_BESA_AVR
swf.npnt  = Npts;
swf.tsb   = TSB;
swf.di    = DI;
swf.sb    = SB;
swf.label = label;
swf.data  = dat;

% FIXME, here combine the channels of the regional sources


