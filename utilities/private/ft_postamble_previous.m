% FT_POSTAMBLE_PREVIOUS

global ft_default

% remember the cfg history of the input data structures
for tmpindx=1:length(ft_default.postamble)
  tmpvar = eval(ft_default.postamble{tmpindx});
  if isfield(tmpvar, 'cfg')
    cfg.previous{tmpindx} = getfield(tmpvar, 'cfg');
  else
    cfg.previous{tmpindx} = [];
  end
end
clear tmpindx tmpvar

if length(cfg.previous)==1
  % replace the cell-array by the single struct
  cfg.previous = cfg.previous{1};
end