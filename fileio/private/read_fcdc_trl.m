function [trl] = read_fcdc_trl(fn);

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
% $Log: read_fcdc_trl.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.2  2008/11/12 16:59:09  roboos
% open explicitely as text using fopen 'rt'
%
% Revision 1.1  2005/04/18 13:47:05  roboos
% added some old and infrequently used functions to the cvs repository
%
% Revision 1.1  2005/04/18 13:43:34  roboos
% included some old functions in the cvs repository, this ensures consistency of the functions between the different network locations
%
% Revision 1.1  2003/04/01 06:54:04  roberto
% *** empty log message ***
%

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

