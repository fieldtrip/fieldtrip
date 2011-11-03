% FT_PREAMBLE_LOADVAR is a helper script that optionally loads one or
% multiple fieldtrip data structures from mat files on disk, as an
% alternative to the user specifying the data structures as input variables
% to the calling function. This makes use of the cfg.inputfile variable.
% 
% Use as
%   ft_preamble loadvar data
%   ft_preamble loadvar source mri
%   ft_preamble loadvar varargin

% use an anonymous function
assign = @(var, val) assignin('caller', var, val);

% the name of the variable is passed in the preamble field
global ft_default

if isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  
  % the input data should be read from file
  if (nargin>1)
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    if isfield(cfg, 'inputlock') && ~isempty(cfg.inputlock)
      mutexlock(cfg.inputlock);
    end
    
    if isequal(ft_default.preamble, {'varargin'}) && ~iscell(cfg.inputfile)
      % this should be a cell-array, oterwise it cannot be assigned to varargin
      cfg.inputfile = {cfg.inputfile};
    end
    
    if iscell(cfg.inputfile)
      if isequal(ft_default.preamble, {'varargin'})
        % read multiple inputs and copy them into varargin
        tmp = {};
        for i=1:length(cfg.inputfile)
          tmp{i} = loadvar(cfg.inputfile{i}, 'data');
        end % for
        assign('varargin', tmp);
        clear i tmp
      else
        % ft_default.preamble is a cell-array containing the variable names
        for i=1:length(cfg.inputfile)
          assign(ft_default.preamble{i}, loadvar(cfg.inputfile{i}, ft_default.preamble{i}));
        end % for
        clear i
      end
    else
      % ft_default.preamble{1} contains the variable name
      assign(ft_default.preamble{1}, loadvar(cfg.inputfile, ft_default.preamble{1}));
    end
    
    if isfield(cfg, 'inputlock') && ~isempty(cfg.inputlock)
      mutexunlock(cfg.inputlock);
    end
  end
end
