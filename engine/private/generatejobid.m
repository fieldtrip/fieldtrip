function job = generatejobid(batch)

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

% remember the current job number for the current batch
jobNum(batch) = job;

