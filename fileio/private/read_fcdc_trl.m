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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

