function [varargout] = funname(varargin)

% READ_CTF_SHM reads metainformation or selected blocks of data from
% shared memory. This function can be used for real-time processing of
% data while it is being acquired.
%
% Use as
%   [msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;
% or
%   [data] = read_ctf_shm(msgNumber);
%   [data] = read_ctf_shm(msgNumber, numValues);
%
% See also WRITE_CTF_SHM

% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);
  mex(mexsrc);
  cd(pwdir);
  success = true;

catch
  % compilation failed
  disp(lasterr);
  error('could not locate MEX file for %s', mexname);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was juist created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end

