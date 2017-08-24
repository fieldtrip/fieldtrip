function id = generatejobid(batch, batchid)

% GENERATEJOBID generates a unique identifier for a job to be submitted to the
% batch queueing system. It maintains an internal counter to allow it to be
% called from multiple qsubfeval instances without the user having to keep
% track of the numbers.
%
% Use as
%   jobid     = generatejobid(batch)
%   batchid   = generatebatchid(batch)
%   sessionid = generatesessionid()
%
% The result is a string like
%   user_host_pid_bM_jN  % as jobid
%   user_host_pid_bM     % as batchid
%   user_host_pid        % as sessionid
% where M is the batch number and N the sequential job number (per batch).
%
% Besides specifying a batch number, it is also possible to specify the full
% batch name as string. This allows the user to override the default
% user_host_pid_bM part of the jobid. This works like
%   jobid     = generatejobid(batch, batchid)
% where for example generatejobid(1, 'freqanalysis') returns the sequence
% of strings 'freqanalysis_j001', 'freqanalysis_j002', ...
%
% See also GENERATEBATCHID, GENERATESESSIONID

% Copyright (C) 2011-2012, Robert Oostenveld
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

% jobNum will be a vector with one entry for each batch number
persistent jobNum

if nargin<1
  error('incorrect number of input arguments');
end

if nargin<2
  batchid = [];
end

if length(jobNum)>=batch
  % increment the previous value
  job = jobNum(batch)+1;
else
  % start with job number one
  job = 1;
end

if ~isempty(batchid)
  id = sprintf('%s_j%03d', batchid, job);
else
  id = sprintf('%s_%s_p%d_b%d_j%03d', getusername(), gethostname(), getpid(), batch, job);
end

% ensure that it can be used as filename, struct fieldname, etc.
id = fixname(id);

% remember the current job number for the current batch
jobNum(batch) = job;
