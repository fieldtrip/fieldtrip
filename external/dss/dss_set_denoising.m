function state = dss_set_denoising(state, denf_function, denf_params)

% Initialize denoising function to state structure
%   [state] = dss_set_denoising(state, denoise_function, denoise_params)
%     state       DSS algorithm state
%     denf_function  Initialized denoising function.
%                    (default:@denoise_fica_tanh)
%     denf_params    Parameters for denoising function. If left out
%                    initializes default parameters for that function.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
  % Initialize default denoising function
  denf_function = @denoise_fica_tanh;
end

if strcmp(class(denf_function),'function_handle')==0
  error('denoise_function parameter must be a function handle');
end

if nargin<3
  % Initialize default parameters for denoising function
  denf_params=struct([]);
  params = feval(denf_function);
  for i=1:length(params.param)
    % for ML65+: denf_params.(params.param{i}) = params.param_value{i};
    eval(strcat('denf_params(1).',params.param{i},'=','params.param_value{i};'));
  end
end

state.denf.h = denf_function;
state.denf.params = denf_params;

% remove alpha and beta
if isfield(state, 'alphaf')
  state = rmfield(state, 'alphaf');
end
if isfield(state, 'betaf')
state = rmfield(state, 'betaf');
end
