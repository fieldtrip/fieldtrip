function retval = defaultbackend(userbackend)

% DEFAULTBACKEND returns a string with the computational backend
% to be used by default. It is determined by looking at various
% environment variables or can be specified by the user.

persistent backend

% check if torque or sge is present and running
if ~isempty(getenv('SGE_ROOT'))
  systembackend = 'sge';
elseif ~isempty(getenv('TORQUEHOME')) || ~isempty(getenv('PBS_VERSION'))
  systembackend = 'torque';
elseif ~isempty(getenv('CONDOR_ARCH'))
  % this has not been tested and I am not 100% sure that this is the right variable to probe
  systembackend = 'condor';
elseif ~isempty(getenv('SLURM_ENABLE')) || ~isempty(getenv('SLURM_SUBMIT_HOST'))
  systembackend = 'slurm';
elseif ~isempty(getenv('LSF_ENVDIR'))
  systembackend = 'lsf';
else
  % backend=local  causes the job to be executed in this MATLAB instance using feval
  % backend=system causes the job to be executed on the same computer using system('matlab -r ...')
  systembackend = 'local';
end

if isempty(backend)
  % set the backend for the first time
  if nargin==0
    backend = systembackend;
  else
    backend = userbackend;
  end
end

if nargin>0 && ~isequal(userbackend, backend)
  warning('switching to the %s backend', userbackend);
  backend = userbackend;
end

% return the backend to the calling function
retval = backend;

if ~mislocked && ~isequal(systembackend, backend)
  mlock
end
