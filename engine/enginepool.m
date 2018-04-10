function retval = enginepool(varargin)

% ENGINEPOOL manages the pool of MATLAB engine workers that is available
% for distributed computing
%
% Use as
%   enginepool open <number> <command>
%   enginepool close
%   enginepool info
%
% The number specifies how many MATLAB engines should be started. In general
% it is advisable to start as many engines as the number of CPU cores.
%
% The command is optional. It can be used to specify the MATLAB version
% and the command-line options. The default for Linux is
%   command = "matlab -singleCompThread -nodesktop -nosplash"
%
% See also ENGINECELLFUN, ENGINEFEVAL, ENGINEGET

% Some undocumented features of this function allow it to maintain an
% internal list that maps between job IDs and engine numbers.
% This use is reserved for engcellfun.
%   enginepool('info')
%   enginepool('block',   index, jobid)
%   enginepool('release', index)
%   enginepool('find',    jobid)

% -----------------------------------------------------------------------
% Copyright (C) 2012, Robert Oostenveld
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

persistent pool

if nargin<1
  cmd = 'info';
else
  cmd = varargin{1};
end

poolsize = length(pool);

switch cmd
  case 'info'
    retval = pool;
    if nargout==0
      isbusy = false(1,numel(pool));
      hasjob = false(1,numel(pool));
      for i=1:numel(pool)
        isbusy(i) = engine('isbusy', i);
        hasjob(i) = ~isempty(pool{i});
      end
      % give some feedback to the user
      if poolsize>0
        fprintf('there are %d engines\n', poolsize);
        fprintf('%d of them have a job\n', sum(hasjob));
        fprintf('%d of them are busy\n',   sum(isbusy));
      else
        fprintf('there are no engines\n');
      end
    end
    
  case 'open'
    if poolsize>1
      error('you first have to close the existing pool before opening another one');
    end
    
    desired = varargin{2};
    if ischar(desired)
      desired = str2double(desired);
    end
    
    if nargin>2
      matlabcmd = varargin{3};
    else
      matlabcmd = [fullfile(matlabroot,'bin','matlab') ' -singleCompThread -nodesktop -nosplash'];
    end
    
    % start the engines
    engine('open', desired, matlabcmd);
    % create a list with job IDs
    pool = cell(1,desired);
    % prevent this function with its persistent variable from getting cleared
    mlock;
    
  case 'close'
    % allow this function with its persistent variable to be cleared
    munlock;
    % clear the list with job IDs
    pool = {};
    % close the engines
    try
      engine('close');
    catch
      % this happens if enginepool and engine get out of sync
      warning(lasterr);
    end
    
  case 'block'
    index = varargin{2};
    if ischar(index)
      index = str2double(index);
    end
    if index>poolsize
      error('invalid index %s, the pool only contains %d workers', index, poolsize);
    end
    
    % add the specified job ID to the persistent list
    jobid = varargin{3};
    pool(index) = {jobid};
    
  case 'release'
    index = varargin{2};
    if ischar(index)
      index = str2double(index);
    end
    if index>poolsize
      error('invalid index %s, the pool only contains %d workers', index, poolsize);
    end
    
    % remove the job ID from the persistent list
    pool(index) = {[]};
    
  case 'find'
    jobid = varargin{2};
    if ischar(jobid)
      jobid = str2double(jobid);
    end
    
    % find the specified job ID in the persistent list
    retval = [];
    for i=1:length(pool)
      if isequal(pool{i}, jobid)
        retval = i;
        break
      end
    end
    
  otherwise
    error('unsupported command "%s"', cmd);
end

if nargout<1
  clear retval
end
