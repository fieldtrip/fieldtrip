function id = generatejobid(batch)

% GENERATEJOBID generates a unique identifier for a job to be submitted to
% the batch queueing system. It consists of:
%
%   user_host_pN_bM_jL
%
% where N,M,L are the calling matlab's process ID, batch number (to be provided
% as an argument) and the sequential job number (per batch), respectively.
%
% generatejobid() without any arguments generates just the first part of a
% possible job ID, so user_host_pN.

persistent jobNum; % jobnum will be numQueues X 1 vector

user = getusername();
host = gethostname();

if (nargin>0 && batch > 0)
  
  if isempty(jobNum)
    jobNum = zeros(batch,1);
    jobNum(batch) = 1;
  elseif numel(jobNum)<batch
    jobNum(batch) = 1;
  else
    jobNum(batch) = jobNum(batch)+1;
  end

  id = sprintf('%s_%s_p%d_b%d_j%03d', user, host, getpid(), batch, jobNum(batch));
  
else
  % batch<0 (or batch=NaN) is used to generate not a real job ID, but just
  % the first part of it (see qsublisten)
  id = [user '_' host '_p' num2str(getpid())];
end
