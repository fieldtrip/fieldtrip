function [trl] = read_fcdc_trl(fn)

% READ_FCDC_TRL reads trial definitions from a file
%
% Given a file which defines N trials, this function returns a Nx3
% matrix with the begin latency, end latency, and the latency offset
% of the first sample of each trial. The latencies are in seconds.
%
% [trl] = read_fcdc_trl(filename)
%
% An FCD trial definition file is formatted like
%   begin   end     offset
%   0.0000  1.0000  0.0000
%   3.0000  4.0000  0.0000
%   5.0000  5.5000  0.0000
%   ... 
% 
% The trial begin and end are given in seconds relative to the start
% of the recorded datafile. The offset is given in seconds and indicates
% the latency of the first sample, relative to the trial marker or
% trigger. E.g., given a trigger at 7000ms (relative to the recording
% begin), a trial of 1000ms with a pretrigger interval of 300ms would
% correspond to "6.700 7.700 -0.300".

% Copyright (C) 2003, Robert Oostenveld
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

fid = fopen(fn, 'rt');
if fid<0
   error('could not open file');
end

trl = [];
while ~feof(fid)
  tmp = fscanf(fid, '%f %f %f', 3);
  if ~isempty(tmp)
    trl = [trl; tmp'];
  end
end

fclose(fid);  

