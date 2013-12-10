function opt = bsearch_opt(opt)
%BSEARCH_OPT Default options for backward search
%
%   Default options for BSEARCH
%
%    opt=bsearch_options;
%      return default otions
%    opt=bsearch_options(opt);
%      fill empty options with default values
%
%   The options and defaults are
%   display (0)
%     0 off
%     1 on 
%     2 more
%   nsel(1)
%     minimum number of columns, Inf means all
%   stop (0)
%     stop if value does not increase
%   stopvalue (-Inf)
%     stop if the current best value smaller than stopvalue

%	Copyright (c) Aki Vehtari (2004-2007)

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


if nargin < 1
  opt=[];
end

if ~isfield(opt,'display')
  opt.display=0;
end
if ~isfield(opt,'nsel')
  opt.nsel=1;
end
if ~isfield(opt,'stop')
  opt.stop=0;
end
if ~isfield(opt,'stopvalue')
  opt.stopvalue=-Inf;
end
