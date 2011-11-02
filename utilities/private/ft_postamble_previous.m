% FT_POSTAMBLE_PREVIOUS adds the cfg structure from the input data
% structure(s) to the current configuration structure as cfg.previous
%
% Use as
%   ft_postamble previous inputvar
%   ft_postamble previous inputvar1 inputvar2
%   ft_postamble previous varargin

global ft_default

% remember the cfg history of the input data structures
cfg.previous = {};
for tmpindx=1:length(ft_default.postamble)
  if exist(ft_default.postamble{tmpindx}, 'var')
    tmpvar = eval(ft_default.postamble{tmpindx});
  else
    tmpvar = [];
  end
  if isfield(tmpvar, 'cfg')
    cfg.previous{tmpindx} = tmpvar.cfg;
  else
    cfg.previous{tmpindx} = [];
  end
end
clear tmpindx tmpvar

if length(cfg.previous)==1
  % replace the cell-array by the single struct
  cfg.previous = cfg.previous{1};
end
