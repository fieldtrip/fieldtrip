function opt = metrop2_opt(opt)
%METROP2_OPT Default options for Metropolis sampling.
%
%   Default options for METROP
%
%    opt=metrop2_opt;
%      return default options
%    opt=metrop2_opt(metropopt);
%      fill empty options with default values
%
%   The options and defaults are
%   display (0)
%     1 to display the energy values and rejection threshold at
%       each step of the Markov chain
%     2 to display also position vectors at each step
%   nsamples (1)
%     the number of samples retained from the Markov chain
%   nomit (0)
%     the number of samples omitted from the start of the chain
%   stddev (0.9)
%     the variance of the proposal distribution; default 1.
%
%	See also
%	METROP2

%	Copyright (c) 1998-2000 Aki Vehtari 

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.
if nargin < 1
  opt=[];
end

if ~isfield(opt,'display')
  opt.display=0;
end
if ~isfield(opt,'nsamples') | opt.nsamples < 1
  opt.nsamples=1;
end
if ~isfield(opt,'nomit') | opt.nomit < 0
  opt.nomit=0;
end
if ~isfield(opt,'stddev') | opt.stddev < 0
  opt.stddev=0.1;
end
