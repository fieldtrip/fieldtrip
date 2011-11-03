% FT_POSTAMBLE_SAVEVAR is a helper script that optionally saves the output
% fieldtrip data structures to a mat file on disk. This is useful for
% batching and for distributed processing. This makes use of the
% cfg.outputfile variable.
%
% Use as
%   ft_preamble savevar data
%   ft_preamble savevar source mri

global ft_default

% the output data should be saved to a MATLAB file
if isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile)
  
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
