function [stat] = statfun_mean(cfg, dat, design)

% STATFUN_MEAN does not depend on the design matrix

stat = mean(dat,2);

