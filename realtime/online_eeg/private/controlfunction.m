function controlfunction(cfg, varargin)

% CONTROLFUNCTION is a helper function that controls external devides
% and that can provide neurofeedback.
%
% Use as
%   controlfunction(cfg, value1, value2, ...)
% where the configuration should contain
%   cfg.method = string
%
% Depending on the control method, additional configuration fields and one
% or multiple input values are required. Please look at the code for all
% details.
%
% See also FT_WRITE_EVENT

% get the method, the default is to display all values
cfg.method = ft_getopt(cfg, 'method', 'disp');

switch cfg.method
  case 'disp'
    % this can be used for an arbitrary number of input values
    fprintf('controlfunction called\n');
    for i=1:length(varargin)
      fprintf('value%d = \n', i);
      if isempty(varargin{i})
        fprintf('[]');
      else
        disp(varargin{i});
      end
    end
    
  case 'midivolume'
    % requires a single value between 0 and 1
    if length(varargin)~=1 || varargin{1}<0 || varargin{1}>1
      error('invalid input');
    end
    
  case 'infrared_tv_on'
    % requires no input values
    if ~isempty(varargin)
      error('invalid input');
    end
    
  case 'infrared_tv_off'
    % requires no input values
    if ~isempty(varargin)
      error('invalid input');
    end
    
  otherwise
    error('unsupported control method "%s"', cfg.method);
end
