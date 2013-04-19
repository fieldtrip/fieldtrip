function opt = sls1mm_opt(opt)
%SLS1MM_OPT Default options for Slice Sampling
%
%  Default options for SLS1MM
%
%    opt = sls1mm_opt;
%      return default options
%    opt = sls1mm_opt(opt);
%      fill empty options with default values
%
%  The options and defaults are (default value)
%    maxiter (50)
%      maximum number of iterations for the shrinkage procedure, in case
%      this is exceeded the other parameters might not be properly assigned
%    mmlimits ([0; 1])
%      absolute limits for the slice (minmax)
%
%  See also
%    SLS1MM

%       Copyright (c) 2003-2004 Toni Auranen
%       Copyright (c) 2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 1
  opt=[];
end

if ~isfield(opt,'maxiter')
  opt.maxiter = 50;
end
if ~isfield(opt,'mmlimits')
  opt.mmlimits = [0; 1];
end
