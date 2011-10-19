% FT_POSTAMBLE_SAVEVAR

global ft_default

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexlock(cfg.outputlock);
  end
  
  % ft_default.postamble{1} contains the name of the variable
  savevar(cfg.outputfile, ft_default.postamble{1}, eval(ft_default.postamble{1}));
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexunlock(cfg.outputlock);
  end
  
  if ~nargout
    % do not return the output variable "ans"
    clear(ft_default.postamble{1});
  end
end
