function id = generatejobid(qnum)
% GENERATEJOBID generates a unique identifier for a job to be submitted to
% a batch queue. It consists of:
%
%   user@host_pN_qM_jobL
%
% where N,M,L are the calling matlab's process ID, queue number (to be provided
% as an argument), and sequential job number (per queue), respectively.

user = getusername();
if strcmp(user,'<unknown>')
  user = 'unknownuser'; % don't want <> in filenames
end

host = gethostname();

persistent jobNum; % jobnum will be numQueues X 1 vector
if isempty(jobNum)
  jobNum = zeros(qnum,1);
  jobNum(qnum) = 1;
elseif numel(jobNum)<qnum
  jobNum(qnum) = 1;
else
  jobNum(qnum) = jobNum(qnum)+1;
end

% matlab does not like an @ in filenames
id = [user '_' host '_p' num2str(getpid()) '_q' num2str(qnum) '_job' num2str(jobNum(qnum))];

function host = gethostname()
  if (ispc())
    host = getenv('ComputerName');
  else
    host = getenv('HOSTNAME');
  end
  
  host = strtok(host, '.'); % dots in filenames are not allowed by matlab
  
  if (isempty(host))
    host = 'unknownhost';
  end
end

end