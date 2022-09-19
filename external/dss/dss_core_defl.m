function state = dss_core_defl(state)
% Deflation DSS algorithm
%   state = dss_core_defl(state)
%     Calculates denoising source separation separately for each component.
%     Begins or continues calculation defined by state structure and
%     returns the result as state structure. Calculation can be interrupted
%     from the keyboard if JAVA is enabled or if testkeypress MEX
%     function exists. Used CPU time is recorded in the state.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: dss_core_defl.m,v 1.25 2005/12/02 12:23:18 jaakkos Exp $

% -- Internal variables
% wdim            Dimension of the whitened data
% interrupt_iteration
% interrupt_componenet
% start_time
% ...

dss_message(state,2,'Extracting components in deflationary DSS\n');

% -- Initialize local variables
start_time = cputime;
interrupt_iteration = 0;
interrupt_component = 0;
% Dimension of the whitened data
if iscell(state.Y)
  wdim = size(state.Y{1},1);
else
  wdim = size(state.Y, 1);
end
sdim = state.sdim;

% -- Determine if continuing old calculation or starting new
if ~isfield(state,'component'); state.component=0; end
if ~isfield(state,'iteration'); state.iteration=0; end
if state.component>0 && state.component<=sdim
  dss_message(state,1,sprintf('Continuing on component %d, iteration %d.\n', state.component, state.iteration));
else
  state.component=1;
  state.iteration=0;
  if ~isfield(state,'W')
    %state.W = zeros(sdim,wdim);
  end
end

state = dss_check_adaptivity(state);

if ~isfield(state,'W'); state.W = []; end
if ~isfield(state,'S'); state.S = []; end

% -- for each component
while state.component <= sdim
  component = state.component;

  if state.iteration==0
    dss_message(state,1,sprintf('Calculate component %d ',component));

    if size(state.W,1)>=component
      dss_message(state,1,'with predefined w\n');
      state.w = state.W(component,:)';
    elseif size(state.S,1)>=component
      dss_message(state,1,'with predefined s\n');
      state.w = state.Y * state.S(component, :)';
    else
      dss_message(state,1,'with random w\n');
      state.w = randn(wdim, 1);
    end

    % -- orthogonalization
    if component > 1
      W = state.W(1:(component-1),:);
    else 
      W = zeros(wdim,1);
    end
    [state.orthof.params, state.w] = feval(state.orthof.h, state.orthof.params, W, state.w);
    
    % estimate
    state.s = state.w' * state.Y;
    state.iteration = 1;
    state.last_iteration = 0;
  else
    dss_message(state,2,'Continue iteration');
  end

  % -- for each iteration
  dss_message(state,3,'Iterating: ');
  while 1
    % ---- calculate new s
    % -- denoising
    [state.denf.params, s_new] = feval(state.denf.h, state.denf.params, ...
				       state.s, state);

    % -- alpha (normalization)
    if state.adapt_alpha
      [state.alphaf.params, state.alpha] = feval(state.alphaf.h, state.alphaf.params, state);
    end
    if state.alpha~=1; s_new = state.alpha * s_new; end

    % -- beta (spectral shift)
    if state.adapt_beta
      [state.betaf.params, state.beta] = feval(state.betaf.h, state.betaf.params, state);
    end
    if state.beta~=0; s_new = s_new + state.beta * state.s; end

    state.s = s_new;

    % ---- calculate new w
    state.w_old = state.w;
    % -- re-estimate projection
    if iscell(state.Y)
      state.w = sum(state.Y * cellctranspose(state.s), 2);
    else
      state.w = state.Y * state.s';
    end
    
    % -- orthogonalization
    if component > 1
      W = state.W(1:(component-1),:);
    else
      W = zeros(wdim,1); 
    end
    [state.orthof.params, state.w] = feval(state.orthof.h, state.orthof.params, W, state.w);

    % check if sign has changed
    signum = sign(state.w' * state.w_old);
    if ~isnan(signum) 
      if signum~=0
	    state.w = signum * state.w;
      end
      % Show progress
      signumstr='-0+';
      dss_message(state,3,signumstr(signum+2));
    end
      
    % -- learning rate adaptation
    if state.adapt_gamma
      [state.gammaf.params, state.gamma] = ...
          feval(state.gammaf.h, state.gammaf.params, state);
      state.w = state.w_old + ...
	  (state.w - state.w_old) * state.gamma;
    
      % -- orthogonalize adapted projection
      if component > 1
        W = state.W(1:(component-1),:);
      else
        W = zeros(wdim,1);
      end
      [state.orthof.params, state.w] = feval(state.orthof.h, state.orthof.params, W, state.w);
    end
    
    
    % -- update new signal estmate
    state.s = state.w' * state.Y;
    
    % -- check stopping criteria
    [state.last_iteration, state.stopf.params] = ...
	feval(state.stopf.h, state.stopf.params, state);

    % -- report iteration
    if isfield(state, 'reportf')
      for f=1:length(state.reportf)
        state.report = feval(state.reportf{f}, state.report, state);
      end
    end
    if state.last_iteration, break, end

    state.iteration = state.iteration + 1;

    % -- Test keypress interrupt
%     if keyboard_interrupt
%       in = input(['\nENTER ''i'' TO STOP NOW OR ''c'' TO STOP AFTER' ...
% 		  ' COMPONENT IS FOUND. ANY OTHER KEY TO CONTINUE: '],'s');
%       if in =='i'
% 	interrupt_iteration = 1;
% 	break;
%       end
%       if in =='c', interrupt_component = 1; end
%     end
      
  end

  % -- store found component
  state.W(state.component, :) = state.w';
  
  if interrupt_iteration
    disp('Interrupted');
    break;
  end

  dss_message(state,3,'\n');

  state.component = state.component + 1;
  state.iteration = 0;

  if interrupt_component
    disp('Interrupted');
    break;
  end

end

state.S = state.W * state.Y;

% -- record the used cpu time
if ~isfield(state, 'cputime'); state.cputime = 0; end
state.cputime = state.cputime + cputime - start_time;
