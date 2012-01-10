  function num = qsublisten(callback, varargin)

% QSUBLISTEN checks whether jobs, submitted by qsubfeval, have been
% completed. Whenever a job returns, it executes the provided callback function
% (should be a function handle), with the job ID as an input argument. Results
% can then be retrieved by calling QSUBGET. If a cell array is provided as
% a the 'filter' option (see below), the second input argument passed to the
% callback function will be an index into this cell array (to facilitate
% checking which job returned in the callback function).
%
% Note that this function is blocking; i.e., it only returns after a
% certain criterion has been met.
%
% Arguments can be supplied with key-value pairs:
%     maxnum      = maximum number of jobs to collect, function will return
%                   after this is reached. Default = Inf; so it is highly
%                   recommended you provide something here, since with
%                   maxnum=Inf the function will never return.
%     filter      = regular expression filter for job IDs to respond to.
%                   The default tests for jobs generated from the current
%                   matlab process. A cell array of strings can be
%                   provided; in that case, exact match is required.
%     sleep       = number of seconds to sleep between checks (default=0)
%
% This function returns the number of jobs that were collected and for
% which the callback function was called.

maxnum = ft_getopt(varargin, 'maxnum', Inf);
filter = ft_getopt(varargin, 'filter', [generatesesionid '.*']);
sleep  = ft_getopt(varargin, 'sleep', 0);

if ischar(filter)
  regexpFilt = 1;
elseif iscellstr(filter)
  regexpFilt = 0;
else
  error('filter should either be a regexp string or cell array of exact-match strings');
end

% keep track of which job IDs we have already recognized and fired the
% callback for
foundJobs = [];

curPwd = getcustompwd();

num = 0;
while (num < maxnum)
  
  files = dir();
  for k = 1:numel(files)
    
    % preliminary filter to get just the qsub-specific output files
    jobid = regexp(files(k).name, '^(.*)\.o.*$', 'tokens');
    
    if ~isempty(jobid) && isempty(findstr(foundJobs, jobid{1}{1}))
      jobid = jobid{1}{1};
      
      % wait until not only the stdout file exists, but also the stderr and
      % _output.mat. If we fire the callback before all three files are
      % present, a subsequent call to qsubget will fail
      outputfile   = fullfile(curPwd, sprintf('%s_output.mat', jobid));
      logerr       = fullfile(curPwd, sprintf('%s.e*', jobid));
      while ~exist(outputfile,'file') || ~isfile(logerr)
        pause(0.01);
      end
    
      if (regexpFilt && ~isempty(regexp(jobid, filter, 'once'))) || ...
          (~regexpFilt && ~isempty(strmatch(jobid, filter, 'exact')))
        
        if (~regexpFilt && nargin(callback)>1)
          % also provide an index into the filter array
          callback(jobid, strmatch(jobid, filter, 'exact'));
        else
          callback(jobid);
        end
        
        num = num+1;
        foundJobs = [foundJobs '|' jobid];
      end
      
    end
  end
  
  if (sleep > 0)
    pause(sleep);
  end
  
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function that detects a file, even with a wildcard in the filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = isfile(name)
tmp = dir(name);
status = length(tmp)==1 && ~tmp.isdir;
end

end