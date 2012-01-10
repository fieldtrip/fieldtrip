function id = generatejobid(batch)

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
% See also GENERATEBATCHID, GENERATESESSIONID

% jobNum will be numQueues X 1 vector
persistent jobNum

if nargin~=1
  error('incorrect number of input arguments');
end

if length(jobNum)>=batch
  % increment the previous value
  job = jobNum(batch)+1;
else
  % start with job number one
  job = 1;
end

id = sprintf('%s_%s_p%d_b%d_j%03d', getusername(), gethostname(), getpid(), batch, job);

% remember the current job number for the current batch
jobNum(batch) = job;

