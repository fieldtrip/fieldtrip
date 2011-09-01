function id = generatejobid(qnum)

% GENERATEJOBID generates a unique identifier for a job to be submitted
% to a batch queue. It consists of:
%
%   user_host_pN_qM_jL
%
% where N,M,L are the calling matlab's process ID, queue number (to
% be provided as an argument), and sequential job number (per queue),
% respectively.

persistent jobNum; % jobnum will be numQueues X 1 vector

host = gethostname();
user = getusername();

if strcmp(user,'<unknown>')
  user = 'unknownuser'; % don't want <> in filenames
end

if isempty(jobNum)
  jobNum = zeros(qnum,1);
  jobNum(qnum) = 1;
elseif numel(jobNum)<qnum
  jobNum(qnum) = 1;
else
  jobNum(qnum) = jobNum(qnum)+1;
end

% matlab does not like an @ in filenames
id = sprintf('%s_%s_p%d_q%d_j%02d', user, host, getpid(), qnum, jobNum(qnum));

