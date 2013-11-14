function [state, B, A] = denss(X_or_state, parameters)
% Main DSS algorithm
%   [state, B, A] = denss(X or state)
%   [state, B, A] = denss(X or state, [params])
%   [state, B, A] = denss(X or state, {param1, 'value1', param2, 'value2', ...})
%     X           Mixed signals
%     state       Existing or new algorithm state structure
%     parameters  Optional parameters
%     B           Unmixing matrix S = B * X
%     A           Mixing matrix X = A * S 
%
%   Main entry point for DSS algorithm. When called with mixed
%   signals (X) creates state structure. Given optional parameters
%   are inserted into new or existing state structure.
%   Calls algorithm core defined by 'algorithm' state/param
%   variable.
%
%   See dss_create_state for description of state and parameter
%   variables.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

% -- Initialize state
if nargin>=2 | ~isstruct(X_or_state)
  % -- user gave X or parameters, initialize state
  if nargin<2; parameters = []; end
  state = dss_create_state(X_or_state, parameters);
else
  state = X_or_state;
end

% -- Contract input dimensions
state = contract_dim(state, 'X');
state = contract_dim(state, 'Y');
state = contract_dim(state, 'S');

% read the version information
VERSION = fopen('VERSION','r');

if VERSION ~= -1
  version = fgets(VERSION);
  version = version(1:end-1);
  dss_message(state,1,sprintf('Calculating DSS (release v%s)\n',version));
else
  dss_message(state,1,'Calculating DSS (unknown release version)\n');
end  
% -- Preprocessing
if ~isfield(state, 'Y')
  state = dss_preprocess(state.X, state);  
else
  % Sphered data is already available, set sphering to identity
  dss_message(state,2,'Sphered data exists, using it.');
  if ~isfield(state, 'V') | ~isfield(state, 'dV')
    state.V = diag(ones(size(state.Y,1),1));
    state.dV = state.V;
    dss_message(state,2,' No sphering matrix given. Assuming identity matrix.\n');
  else
    % using given sphering matrix
    dss_message(state,2,'\n');
  end
end

state.wdim = size(state.Y, 1);
if ~isfield(state, 'sdim'); state.sdim = state.wdim; end

% -- Call correct DSS algorithm implementation
state = feval(['dss_core_' state.algorithm], state);

% -- return the results
state.B = state.W * state.V;
B = state.B;
state.A = state.dV * state.W';
A = state.A;

%----------------------
function state = contract_dim(state, field_name)
if isfield(state, field_name)
  % for ML65+: dims = size(state.(field_name));
  dims = size(getfield(state, field_name));
  if length(dims)>2
    if isfield(state, 'input_dims') if state.input_dims ~= dims
        error('Input signal dimension mismatch');
    end; end;
    state.input_dims = dims(2:length(dims));
    % for ML65+: state.(field_name) = reshape(state.(field_name), dims(1), []);
    state = setfield(state, field_name, reshape(getfield(state, field_name), dims(1), []));
    dss_message(state,1,sprintf('Contracting [%s] dimensional input %s to [%s] dimension.\n', tostring(dims), field_name, tostring(size(getfield(state,field_name)))));
  end
end

