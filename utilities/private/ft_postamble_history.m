% FT_POSTAMBLE_HISTORY

global ft_default
varname = ft_default.postamble{1};

% remember the configuration detaild in the output variable
eval(sprintf('%s.cfg = cfg;', varname));
