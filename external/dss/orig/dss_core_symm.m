function state = dss_core_symm(state)
% Symmetric DSS algorithm
%   state = dss_core_symm(state)
%     Calculates denoising source separation for a collection of
%     mixed signals.
%     Components are forced orthogonal with symmetric orthogonalization.
%     Begins or continues calculation defined by state structure and
%     returns the result as state structure. Calculation can be interrupted
%     from the keyboard if JAVA is enabled or if testkeypress MEX
%     function exists. Used CPU time is recorded in the state.


% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

dss_message(state,2,'Extracting components in symmetric DSS\n');

% -- Initialize local variables
start_time = cputime;
user_interrupt = 0;
wdim = size(state.Y, 1);
sdim = state.sdim;

if isfield(state,'iteration')
  dss_message(state,1,sprintf('Continuing iteration %d.\n', state.iteration));
else
  state.iteration=1;
  dss_message(state,1,sprintf(...
    'Mixed signal dim: %d, original signal dim: %d\n', wdim, sdim));

  if isfield(state,'W')
    dss_message(state,1,sprintf('Using %d predefined w vectors\n',size(state.W,1)));
  end
  if isfield(state,'S') if size(state.S,1)>size(state.W,1)
    b = size(state.W,1);
    e = size(state.S,1);
    dss_message(state,1,sprintf('Using %d predefined s vectors\n',e-b));
    state.W(b+1:e,:) = state.S(b+1:e,:) * state.Y';
  end; end
  if ~isfield(state, 'W') state.W = []; end;
  state.W(size(state.W,1)+1:sdim,:) = orth(randn(sdim-size(state.W,1),wdim)')';

  % -- orthogonalization
  [state.orthof.params, state.W] = feval(state.orthof.h, state.orthof.params, state.W);
  
  % estimate
  state.S = state.W * state.Y;
end

state = dss_check_adaptivity(state);

% -- for each iteration
while 1
  % ---- calculate new s
  % -- denoising
  [state.denf.params, S_new] = feval(state.denf.h, state.denf.params, ...
				     state.S, state);

  % -- alpha (normalization)
  if state.adapt_alpha
    state.alpha = feval(state.alphaf.h, state.alphaf.params, state);
  end
  if state.alpha~=1
    S_new = repmat(state.alpha, 1, size(S_new, 2)) * S_new;
  end
  
  % -- beta (spectral shift)
  if state.adapt_beta
    state.beta = feval(state.betaf.h, state.betaf.params, state);
  end
  if state.beta~=0
    if length(state.beta)==1
      S_new = S_new + state.beta .* state.S;
    else
      S_new = S_new + repmat(state.beta, 1, size(state.S, 2)) ...
	      .* state.S;
    end
  end
 
  state.S = S_new;
  
  % ---- calculate new w
  state.W_old = state.W;
  % -- re-estimate projection
  state.W = state.S * state.Y';
  
  % -- orthogonalization
  [state.orthof.params, state.W] = feval(state.orthof.h, state.orthof.params, state.W);
  
  % check if sign has changed
  %%% TODO: all the signs should be changed individually
  signum = sign(sum(state.W .* state.W_old, 2)) * ones(1, wdim);
  state.W = signum .* state.W;

  % -- learning rate adaptation
  if state.adapt_gamma
    [state.gammaf.params, gamma] = feval(state.gammaf.h, ...
  					 state.gammaf.params, ...
					 state);
    state.W = state.W_old + (state.W - state.W_old) ...
	.* repmat(gamma, size(state.W, 1), 1);
    % -- orthogonalize adapted projection
    [state.orthof.params, state.W] = feval(state.orthof.h, state.orthof.params, state.W);
  end
  
  % -- update new signal estmate
  state.S = state.W * state.Y;

  % -- check stopping criteria
  [state.last_iteration, state.stopf.params] = ...
      feval(state.stopf.h, state.stopf.params, state);
  
  % -- report iteration
  if isfield(state, 'reportf')
    for f=1:length(state.reportf)
      state.report = feval(state.reportf{f}, state.report, state);
    end
  end

  if state.last_iteration, return, end
  
  state.iteration = state.iteration + 1;
  
  % -- Test keypress interrupt
  if keyboard_interrupt
    in = input('\nENTER ''i'' TO STOP NOW. ANY OTHER KEY TO CONTINUE: ','s');
    if in =='i'
      user_interrupt = 1;
      return;
    end
  end

  dss_message(state,3,'.');

end

if user_interrupt
  disp('Interrupted');
end

dss_message(state,3,'\n');

% -- record the used cpu time
if ~isfield(state, 'cputime'); state.cputime = 0; end
state.cputime = state.cputime + cputime - start_time;

