function [meg] = read_ctf_dat(filename)

% READ_CTF_DAT reads MEG data from an ascii format CTF file
%
% meg = read_ctf_dat(filename)
% 
% returns a structure with the following fields:
%   meg.data        Nchans x Ntime
%   meg.time        1xNtime in miliseconds
%   meg.trigger     1xNtime with trigger values
%   meg.label       1xNchans cell array with channel labels (string)

% Copyright (C) 2002, Robert Oostenveld
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

fid = fopen(filename, 'r');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

% read the sample number
line = fgetl(fid);
[tok, rem] = strtok(line, ':');
meg.sample = str2num(rem(2:end));

% read the time of each sample and convert to miliseconds
line = fgetl(fid);
[tok, rem] = strtok(line, ':');
meg.time = 1000*str2num(rem(2:end));

% read the trigger channel 
line = fgetl(fid);
[tok, rem] = strtok(line, ':');
meg.trigger = str2num(rem(2:end));

% read the rest of the data
meg.data = [];
meg.label = {};
chan = 0;
while (1)
  line = fgetl(fid);
  if ~isempty(line) && line==-1
    % reached end of file
    break
  end
  [tok, rem] = strtok(line, ':');
  if ~isempty(rem)
    chan = chan + 1;
    meg.data(chan, :) = str2num(rem(2:end));
    meg.label{chan} = fliplr(deblank(fliplr(deblank(tok))));
  end
end

% convert to fT (?)
meg.data = meg.data * 1e15;

% apparently multiple copies of the data can be stored in the file
% representing the average, the standard deviation, etc.
% keep only the first part of the data -> average
tmp = find(diff(meg.time)<0);
if ~isempty(tmp)
  meg.data = meg.data(:,1:tmp(1));
  meg.time = meg.time(1:tmp(1));
  meg.trigger = meg.trigger(1:tmp(1));
  meg.sample = meg.sample(1:tmp(1));
end

