function backend = defaultbackend

% DEFAULTBACKEND returns a string with the computational backend
% to be used by default. It is determined by looking at various
% environment variables.

% check if torque or sge is present and running
if ~isempty(getenv('SGE_ROOT'))
  backend = 'sge';
elseif ~isempty(getenv('TORQUEHOME')) || ~isempty(getenv('PBS_VERSION'))
  backend = 'torque';
elseif ~isempty(getenv('CONDOR_ARCH'))
  % this has not been tested and I am not 100% sure that this is the right variable to probe
  backend = 'condor';
elseif ~isempty(getenv('SLURM_ENABLE'))
  backend = 'slurm';
elseif ~isempty(getenv('LSF_ENVDIR'))
  backend = 'lsf';
else
  % backend=local causes the job to be executed in this MATLAB by feval
  % backend=system causes the job to be executed on the same computer using system('matlab -r ...')
  backend = 'local';
end


