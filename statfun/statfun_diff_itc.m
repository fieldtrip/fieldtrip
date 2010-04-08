function [s] = statfun_xxx(cfg, dat, design);

% STATFUN_xxx is a function for computing a statistic for the relation
% between biological data and a design vector containing trial
% classifications or another independent variable
%
% This function is called by STATISTICS_RANDOM, where you can specify
% cfg.statistic = 'xxx' which will be evaluated as statfun_xxx.
%
% The external interface of this function has to be
%   [s] = statfun_xxx(cfg, dat, design);
% where
%   dat    contains the biological data, Nvoxels x Nreplications
%   design contains the independent variable,  1 x Nreplications
%
% Additional settings can be passed through to this function using
% the cfg structure.
%   cfg.complex       = 'absdiff', 'diffabs'

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the defaults
if ~isfield(cfg, 'complex'), cfg.complex = 'diffabs';   end

selA = find(design(cfg.ivar,:)==1);
selB = find(design(cfg.ivar,:)==2);
dfA  = length(selA);
dfB  = length(selB);
if (dfA+dfB)<size(design, 2)
  error('inappropriate design, should contain 1''s and 2''s');
end
if isreal(dat)
  error('the data should be complex, i.e. computed with freqanalysis and cfg.output="fourier"');
end
% normalise the complex data in each trial
dat = dat./abs(dat);

if strcmp(cfg.complex, 'diffabs')
  % first compute the absolute, then take teh difference
  % this is not sensitive to phase differences
  itcA = abs(mean(dat(:,selA), 2)); % ITC is the length of the average complex numbers
  itcB = abs(mean(dat(:,selB), 2)); % ITC is the length of the average complex numbers
  s = itcA - itcB;
else
  % first compute the difference, then take teh absolute
  % this is sensitive to phase differences
  itcA = mean(dat(:,selA), 2); % ITC is here the average complex number
  itcB = mean(dat(:,selB), 2); % ITC is here the average complex number
  s = abs(itcA - itcB);
end

s.stat = s;

