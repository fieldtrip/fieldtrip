function [argout, optout] = qsubexec(jobid)

% QSUBEXEC is a helper function to execute a job on the Torque, SGE, PBS
% or SLURM batch queue system. Normally you should not start this function
% yourself, but rather use QSUBCELLFUN or QSUBFEVAL.
%
% This function performs the following tasks
% - load the function name, input arguments and further options from the input file
% - evaluate the desired function on the input arguments using PEEREXEC
% - save the output arguments to an output file
%
% This function should be started from the linux command line as follows
%   qsub /opt/bin/matlab -r "qsubexec(jobid); exit"
% which starts the MATLAB executable, executes this function and exits
% MATLAB to leave your batch job in a clean state. The jobid is
% automatically translated into the input and output file names, which
% have to reside on a shared network file system.
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

try
  [p, jobid] = fileparts(jobid);
  if isempty(p)
    p = pwd();
  end

  inputfile  = fullfile(p, sprintf('%s_input.mat',   jobid));
  outputfile = fullfile(p, sprintf('%s_output.mat_', jobid));

  stopwatch = tic;
  while (~exist(inputfile, 'file') && toc(stopwatch)<60)
    % the underlying NFS file system might be slow in updating, wait for up to 60 seconds
    warning('the input file %s does not yet exist', inputfile);
    pausejava(10);
  end
  clear stopwatch
  
  if ~exist(inputfile, 'file')
    error('timeout while waiting for the input file %s', inputfile);
  end
  
  % the input file contains a function handle
  % catch the warning if the function handle is not recognized
  lastwarn('');
  tmp = load(inputfile);
  [lastmsg, lastid] = lastwarn;
  
  if isequal(lastid, 'MATLAB:dispatcher:UnresolvedFunctionHandle')
    % argin{1} or argin{2} might be a private function
    whichfunction = ft_getopt(tmp.optin, 'whichfunction');
    if ~isempty(whichfunction) && exist(whichfunction, 'file')
      warning('assuming %s as full function name', whichfunction);
      oldpwd = pwd;
      [fundir, funname] = fileparts(whichfunction);
      cd(fundir)
      if isequal(tmp.argin{1}, @cellfun) || isequal(tmp.argin{1}, 'cellfun')
        tmp.argin{2} = str2func(funname);
      else
        tmp.argin{1} = str2func(funname);
      end
      cd(oldpwd);
    end
  end
  
  rerunable = ft_getopt(tmp.optin, 'rerunable', 'no');
  rerunable = istrue(rerunable);

  if ~rerunable
    % delete it as soon as possible to avoid the working directory from getting cluttered
    delete(inputfile);
  end
  
  argin = tmp.argin; % this includes the function name and the input arguments
  optin = tmp.optin; % this includes the path setting, the pwd, the global variables, etc.
  
  [argout, optout] = fexec(argin, optin);
  
  % if variables < ~500 MB, store it in old (uncompressed) format, which is faster
  s1 = whos('argout');
  s2 = whos('optout');
  if (s1.bytes + s2.bytes < 500000000)
    save(outputfile, 'argout', 'optout', '-v6');
  else
    save(outputfile, 'argout', 'optout', '-v7.3');
  end

  if rerunable
    % delete it only once the result has been written to disk
    delete(inputfile);
  end
  
  % remove the _ at the end, note that the rename command here is a private mex file
  retval = rename(outputfile, outputfile(1:end-1));
  if retval~=0
    error('problem renaming output file %s', outputfile);
  end
  
catch err
  % this is to avoid MATLAB from hanging in case fexec fails, since
  % after the job execution we want MATLAB to exit
  disp(err);
  warning('an error was caught');
  
end % try-catch
