% FT_PREAMBLE_LOADVAR

% use an anonymous function
assign = @(var, val) assignin('caller', var, val);

% the name of the variable is passed in the preamble field
global ft_default

% load optional given inputfile as data
hasdata = (nargin>1);

if isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    if isfield(cfg, 'inputlock') && ~isempty(cfg.inputlock)
      mutexlock(cfg.inputlock);
    end
    
    % ft_default.preamble{1} contains the variable name
    assign(ft_default.preamble{1}, loadvar(cfg.inputfile, ft_default.preamble{1}));
    
    if isfield(cfg, 'inputlock') && ~isempty(cfg.inputlock)
      mutexunlock(cfg.inputlock);
    end
  end
end
