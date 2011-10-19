% FT_POSTAMBLE_SAVEVAR

global ft_default
varname = ft_default.postamble{1};

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, varname, eval(varname));
  if ~nargout
    % do not return the output variable "ans"
    clear(varname);
  end
end

clear varname
