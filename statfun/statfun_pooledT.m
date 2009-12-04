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

selA = find(design(cfg.ivar,:)==1);
selB = find(design(cfg.ivar,:)==2);
dfA  = length(selA);
dfB  = length(selB);
if (dfA+dfB)<size(design, 2)
  error('inappropriate design, should contain 1''s and 2''s');
end
sumA = sum(dat(:,selA), 2);
sumB = sum(dat(:,selB), 2);
s = (sumA - sumB)./sqrt(dfA);

s.stat = s;

