function retval = qsublist(cmd, jobid, pbsid)

% QSUBLIST is a helper function that is used to keep track of all the jobs in a
% submitted batch. specifically, it is used to maintain the mapping between the
% job identifier in the batch queueing system and MATLAB.
%
% Use as
%   qsublist('list')
%   qsublist('killall')
%   qsublist('kill', jobid)
%   qsublist('getjobid', pbsid)
%   qsublist('getpbsid', jobid)
%
% The jobid is the identifier that is used within MATLAB for the file names,
% for example 'roboos_mentat242_p4376_b2_j453'.
%
% The pbsid is the identifier that is used within the batch queueing system,
% for example '15260.torque'.
%
% The following commands can be used by the end-user.
%   'list'      display all jobs
%   'kill'      kill a specific job, based on the jobid
%   'killall'   kill all jobs
%   'getjobid'  return the mathing jobid, given the pbsid
%   'getpbsid'  return the mathing pbsid, given the jobid
%
% The following low-level commands are used by QSUBFEVAL and QSUBGET for job
% maintenance and monitoring.
%   'add'
%   'del'
%   'completed'
%
% See also QSUBCELLFUN, QSUBFEVAL, QSUBGET

% -----------------------------------------------------------------------
% Copyright (C) 2011-2016, Robert Oostenveld
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/
%
% $Id$
% -----------------------------------------------------------------------

persistent list_jobid list_pbsid

% this function should stay in memory to keep the persistent variables for a long time
% locking it ensures that it does not accidentally get cleared if the m-file on disk gets updated
mlock

if ~isempty(list_jobid) && isequal(list_jobid, list_pbsid)
  % it might also be system, but torque, sge, slurm and lsf will have other job identifiers
  backend = 'local';
else
  % use the environment variables to determine the backend
  backend = defaultbackend;
end

if nargin<1
  cmd = 'list';
end

if nargin<2
  jobid = [];
end

if nargin<3
  pbsid = [];
end

if isempty(jobid) && ~isempty(pbsid)
  % get it from the persistent list
  sel = find(strcmp(pbsid, list_pbsid));
  if length(sel)==1
    jobid = list_jobid{sel};
  else
    warning('cannot determine the jobid that corresponds to pbsid %s', pbsid);
  end
end

if isempty(pbsid) && ~isempty(jobid)
  % get it from the persistent list
  sel = find(strcmp(jobid, list_jobid));
  if length(sel)==1
    pbsid = list_pbsid{sel};
  else
    warning('cannot determine the pbsid that corresponds to jobid %s', jobid);
  end
end

switch cmd
  case 'add'
    % add it to the persistent lists
    list_jobid{end+1} = jobid;
    list_pbsid{end+1} = pbsid;

  case 'del'
    sel = strcmp(jobid, list_jobid);
    % remove the job from the persistent lists
    list_jobid(sel) = [];
    list_pbsid(sel) = [];

  case 'kill'
    sel = strcmp(jobid, list_jobid);
    if any(sel)
      % remove it from the batch queue
      switch backend
        case 'torque'
          system(sprintf('qdel %s', pbsid));
        case 'sge'
          system(sprintf('qdel %s', pbsid));
        case 'slurm'
          system(sprintf('scancel --name %s', jobid));
        case 'lsf'
          system(sprintf('bkill %s', pbsid));
        case 'local'
          % cleaning up of local jobs is not supported
        case 'system'
          % cleaning up of system jobs is not supported
      end
      % remove the corresponing files from the shared storage
      system(sprintf('rm -f %s*', jobid));
      % remove it from the persistent lists
      list_jobid(sel) = [];
      list_pbsid(sel) = [];
    end

  case 'killall'
    if ~isempty(list_jobid)
      % give an explicit warning, because chances are that the user will see messages from qdel
      % about jobs that have just completed and hence cannot be deleted any more
      fprintf('cleaning up all scheduled and running jobs, don''t worry if you see warnings from "qdel"\n');
    end
    % start at the end, work towards the begin of the list
    for i=length(list_jobid):-1:1
      qsublist('kill', list_jobid{i}, list_pbsid{i});
    end

  case 'completed'
    % cmd = 'completed' returns whether the job is completed as a boolean
    %
    % It first determines whether the output files exist. If so, it might be that the
    % batch queueing system is still writing to them, hence the next system-specific
    % check also polls the status of the job. First checking the files and then the
    % job status ensures that we don't saturate the torque server with job-status
    % requests.

    curPwd     = getcustompwd();
    outputfile = fullfile(curPwd, sprintf('%s_output.mat', jobid)); % if the job is aborted to a resource violation, there will not be an output file
    logout     = fullfile(curPwd, sprintf('%s.o*', jobid)); % note the wildcard in the file name
    logerr     = fullfile(curPwd, sprintf('%s.e*', jobid)); % note the wildcard in the file name

    % poll the job status to confirm that the job truely completed
    if isfile(logout) && isfile(logerr) && ~isempty(pbsid)
      % only perform the more expensive check once the log files exist
      switch backend
        case 'torque'
          [dum, jobstatus] = system(['qstat ' pbsid ' -f1 | grep job_state | grep -o "= [A-Z]" | grep -o "[A-Z]"']);
          if isempty(jobstatus)
            warning('cannot determine the status for pbsid %s', pbsid);
            retval = 1;
          else
            retval = strcmp(strtrim(jobstatus) ,'C');
          end
        case 'lsf'
          [dum, jobstatus] = system(['bjobs ' pbsid ' | awk ''NR==2'' | awk ''{print $3}'' ']);
          retval = strcmp(strtrim(jobstatus), 'DONE');
        case 'sge'
          [dum, jobstatus] = system(['qstat -s z | grep ' pbsid ' | awk ''{print $5}''']);
          retval = strcmp(strtrim(jobstatus), 'z') | strcmp(strtrim(jobstatus), 'qw');
        case 'slurm'
          % only return the status based on the presence of the output files
          % FIXME it would be good to implement a proper check for slurm as well
          retval = 1;
        case {'local','system'}
          % only return the status based on the presence of the output files
          % there is no way polling the batch execution system
          retval = 1;
      end
    elseif isfile(logout) && isfile(logerr) && isempty(pbsid)
      % we cannot locate the job in the PBS/torque backend (weird, but it happens), hence we have to rely on the e and o files
      % note that the mat file still might be missing, e.g. when the job was killed due to a resource violation
      retval = 1;
    else
      retval = 0;
    end

  case 'list'
    for i=1:length(list_jobid)
      fprintf('%s %s\n', list_jobid{i}, list_pbsid{i});
    end

  case 'getjobid'
    % return the mathing jobid, given the pbsid
    retval = jobid;

  case 'getpbsid'
    % return the mathing pbsid, given the jobid
    retval = pbsid;

  otherwise
    error('unsupported command (%s)', cmd);
end % switch

if length(list_jobid)~=length(list_pbsid)
  error('jobid and pbsid lists are inconsistent');
end

if mislocked && isempty(list_jobid) && isempty(list_pbsid)
  % it is now safe to unload the function and persistent variables from memory
  munlock
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function that detects a file, even with a wildcard in the filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = isfile(name)
tmp = dir(name);
status = length(tmp)==1 && ~tmp.isdir;
