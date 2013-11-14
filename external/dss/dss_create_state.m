function state = dss_create_state(X_or_state, parameters)
% Initializes DSS state structure based on existing state
% and given parameters.
%   state = dss_create_state(X or state, parameters)
%     X           Mixed source, matrix with row vectors as signals.
%     state       New or updated persistent algorithm state structure
%     parameters  Optional parameters that are included in the state.
%                 Overwrites parameters in the state variable.
%
% Description of parameter and state variables:
%   algorithm       DSS algorithm
%                     'defl' = Deflation approach (default)
%		      'symm' = Symmetric approach 
%		      'pca'  = Linear PCA DSS
%   sdim            Number of components to be extracted
%   wdim            Dimension of the data after whitening (preprocessing)
%   verbose         Verbosity
%                     0 = no print output
%                     1 = terse print output (default)
%                     2 = more detailed print output 
%                     3 = verbose print output (shows iteration
%                     progress)
%
%   preprocf.h      Preprocessing function (default: pre_sphere)
%   preprocf.params Modifiable parameter structure
%   orthof.h        Orthogonalization function (default: ortho_default)
%   orthof.params   Modifiable parameter structure
%   denf.h          Denoising function (default: denoise_fica_tanh)
%   denf.params     Modifiable parameter structure
%   stopf.h         Stopping criteria function (default: default_stop)
%   stopf.params    Modifiable parameter structure
%
%   reportf         Vector or list of iteration reporting functions.
%   report          Reporting data used by reporting functions
%
%   alphaf.h        Denoise normalization function (default:none)
%   alphaf.params   Modifiable parameter structure
%   betaf.h         Denoise spectral shift function (default: none)
%   betaf.params    Modifiable parameter structure
%   gammaf.h        Adaptive learning rate function (default:none)
%   gammaf.params   Modifiable parameter structure
%
%   alpha           Constant value for alpha (not adaptive, default: 1)
%   beta            Constant value for beta (not adaptive, default: 0)
%   gamma           Constant value for gamma (not adaptive, default: 1)
%
% Internal data and state variables
%   X          Original mixed source
%   Y          Sphered mixed source
%   A          Mixing matrix
%   B          Unmixing matrix (for X)
%   W          Matrix of all/found projection vectors (for Y)
%   s          Estimate of the current source component (deflation)
%   w          Estimate of the current projection vector (deflation)
%   w_old      Projection vector of the previous iteration (deflation)
%   S          Estimate of all the signal components (symmetric)
%   W_old      Projection vector matrix of the previous iteration (symmetric)
%   iteration  Current iteration number
%   component  Index of currently calculated component (deflation)
%   V          Sphering matrix
%   dV         De-sphering matrix
%
% There are three ways to give values for alpha, beta and gamma:
%  1) Constant value (alpha, beta, gamma),
%  2) function handle (alphaf, betaf, gammaf) with
%    a) constant value (adapt_... = 0) or
%    b) adaptive value (adapt_... = 1)

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<1
  error('DSS initialization requires state structure with input data.');
end

% -- Determine calling type
if isstruct(X_or_state)
  state = X_or_state;
else
  state.X = X_or_state;
end

% -- default functions
def_denoise=struct('h', @denoise_fica_tanh,'params', '');
def_stop   =struct('h', @default_stop, ...
    'params', struct('epsilon', 1e-1,'maxiters',200));
def_preproc=struct('h', @pre_sphere);
def_ortho  =struct('h', @ortho_default,  'params', '');
def_report =struct('data', {{}});

% -- parameter & state variables
% Strcture definition: name, description, [default value], [allowed values]
param_def = {
  {'verbose',      'Output verbosity', 1, {0,1,2,3}},
  {'algorithm',    'DSS algorithm type', 'defl', {'defl', 'symm'}},
  {'preprocf',     'Preprocessing function', def_preproc},
  {'orthof',       'Orthogonalization function', def_ortho},
  {'denf',         'Denoising function', def_denoise},
  {'stopf',        'Stoppping criterion function', def_stop},
  {'reportf',      'Vector of iteration reporting functions'},
  {'sdim',         'Source signal dimension'},
  {'wdim',         'Projection dimension'},
  {'report',       'Reporting data', def_report},

  {'alphaf',      'Denoise normalization function'},
  {'betaf',       'Denoise spectral shift function'},
  {'gammaf',      'Adaptive gain function'},
  {'alpha',       'Fixed value for alpha if alphaf is not defined,',1},
  {'beta',        'Fixed value for beta if betaf is not defined,',0},
  {'gamma',       'Fixed value for gamma if gammaf is not defined,',1},
  {'W',           'Initial guess for projection matrix'}
  {'S',           'Initial guess for components'}
  {'Y',           'Pre-sphered data'}
};

if nargin<2
  params = struct([]);
else
  if length(parameters)==0
    params = struct([]);
  else
    if isstruct(parameters)
      dss_message(state,1,'Parse parameter struct\n');
      params = parameters;
    else
      dss_message(state,1,'Parse parameter list\n');
      if mod(size(parameters,2),2)~=0
        error('Parameter list must contain parameter and value pairs.');
      end
    
      for a=1:2:size(parameters,2)
        parameter = parameters{a};
        if ~ischar(parameter)
          error(['Invalid parameter name at ' num2str(a+1)]);
        end
        value = parameters{a+1};
        % for ML65+: params.(parameter) = value;
        params = setfield(params, parameter, value);
      end
    end
  end
end

% -- Initialize given parameters to state variable
for p=1:size(param_def,1)
  paramname = param_def{p}{1};
  description = param_def{p}{2};
  default_value=[];
  if (size(param_def{p},2)>=3); default_value = param_def{p}{3}; end
  allowed_values=[];
  if (size(param_def{p},2)>=4); allowed_values = param_def{p}{4}; end
  
  % -- copy given parameter value to state
  if isfield(params, paramname)
    % for ML65+: state.(paramname)=params.(paramname);
    state = setfield(state, paramname, getfield(params, paramname));
    %dss_message(state,1,['Parameter ' paramname '\n']);
    % Clear parameter from struct so that invalid parameters are detected
    params = rmfield(params, paramname);
  end

  % -- default value
  if ~isfield(state, paramname) & ~isempty(default_value)
    % for ML65+: state.(paramname) = default_value;
    state = setfield(state, paramname, default_value);
  end

  % TODO: check allowed values
  
  % -- check function parameters
  if isfield(state, paramname)
    if isfield(getfield(state, paramname), 'h') & ~isfield(getfield(state, paramname), 'params')
      %dss_message(state,1, sprintf('Creating empty parameter structure for parameter %s.\n', paramname));

      % for ML65+: state.(paramname).params=[];
      p = getfield(state, paramname);
      p.params = [];
      state = setfield(state, paramname, p);
      
    end
  end

  % -- Print info about variable
  if ~isfield(state, paramname) str = '-';
  else str = tostring(getfield(state, paramname));
  end
  dss_message(state, 2, ['  Parameter ' paramname ': ' str '\n']);
  
end

% -- Check invalid parameters
% TODO: should invalid parameters throw an error?

% for ML65+: fields  = fieldnames(params);
if isa(params, 'struct') fields  = fieldnames(params);
else fields = [];
end

if size(fields,1)>0
  for p = size(fields, 1)
    dss_message(state,1,['IGNORING INVALID PARAMETER: ''' fields{p} '''\n']);
  end
end

% other state initializations
if ~isfield(state, 'report')
  state.report = [];
end

