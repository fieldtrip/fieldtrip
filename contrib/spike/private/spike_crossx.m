function [varargout] = spike_crossx(varargin)

% FT_SPIKE_SUB_CROSSX computes cross-correlations between two
% point-processes. 
%
% Use as [C,bins] = ft_spike_sub_crossx(t1,t2,binsize,nbins)
% 
% Inputs:
%   t1 and t2 are two sorted vectors containing spike-times
%   binsize: the binsize for the cross-correlation histogram in seconds
%   nbins: the number of bins to use for the cross-correlation histogram
%
% Ouputs:
%   C is the cross-correlation histogram
%   bins: vector with the times corresponding to the bin centers
%
% Note that this function is implemented as a mex file. If the mex file is 
% missing on your platform, it will make an attempt to automatically 
% compile it
 
% Copyright (C) 2011-2012 Martin Vinck
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  cd(mexdir);
  mex(mexsrc);
  cd(pwdir);
  success = true;
catch
  % compilation failed
  disp(lasterr);
  if varargin{5}==1
    warning('could not MEX file for %s', mexname);
  end
  cd(pwdir);
  success = false;
end

if success 
  % execute the mex file that was just created
  funname = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end
 
