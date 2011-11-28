function [varargout] = spike_crossx(varargin)

% FT_SPIKE_SUB_CROSSX computes cross-correlations between two
% point-processes. 
%
% Use as [C,bins] = ft_spike_sub_crossx(t1,t2,binsize,bins)
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
 
% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc = [funname '.c'];
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
funname = mfilename;
funhandle = str2func(funname);
[varargout{1:nargout}] = funhandle(varargin{:});
end