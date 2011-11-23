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
cnt = 0;
for tmpindx=1:length(ft_default.postamble)
  if exist(ft_default.postamble{tmpindx}, 'var')
    tmpvar = eval(ft_default.postamble{tmpindx});
  else
    tmpvar = [];
  end
  if isa(tmpvar, 'struct')
    % the variable is a data structure
    cnt=cnt+1;
    if isfield(tmpvar, 'cfg')
      cfg.previous{cnt} = tmpvar.cfg;
    else
      cfg.previous{cnt} = [];
    end
  elseif isa(tmpvar, 'cell')
    % the variable is a cell-array (i.e. most likely a varargin)
    for cellindx=1:numel(tmpvar)
      cnt=cnt+1;
      if isfield(tmpvar{cellindx}, 'cfg')
        cfg.previous{cnt} = tmpvar{cellindx}.cfg;
      else
        cfg.previous{cnt} = [];
      end
    end
  end
end
clear tmpindx tmpvar cellindx cnt

if length(cfg.previous)==1
  % replace the cell-array by the single struct
  cfg.previous = cfg.previous{1};
end
