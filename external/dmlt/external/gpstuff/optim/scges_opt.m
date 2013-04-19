function opt = scges_opt(opt)
%SCGES_OPT Default options for scaled conjugate gradient optimization
%
%  Description
%    OPT=SCGES_OPT returns default options for SCGES
%    OPT=SCGES_OPT(OPT); fills empty options with default values
%
%    The options and defaults are
%    display (0)
%      0 to not display anything
%      1 to display just diagnostic messages
%      2 positive integer to show also the function values and
%        validation values every iteration
%    checkgrad (0)
%      1 to check the user defined gradient function
%    maxiter (1000)
%      maximum number of iterations
%    tolfun (1e-6)
%      termination tolerance on the function value
%    tolx (1e-6)
%      termination tolerance on X
%
%  See also
%    SCGES

%	Copyright (c) Aki Vehtari (1998-2010)

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 1
  opt=[];
end

if ~isfield(opt,'display')
  opt.display=0;
end
if ~isfield(opt,'checkgrad')
  opt.checkgrad=0;
end
if ~isfield(opt,'maxiter') | opt.maxiter < 1
  opt.maxiter=100;
end
if ~isfield(opt,'tolfun') | opt.tolfun < 0
  opt.tolfun=1e-6;
end
if ~isfield(opt,'tolx') | opt.tolx < 0
  opt.tolx=1e-6;
end
if ~isfield(opt,'maxfail') 
  opt.maxfail=20;
end
