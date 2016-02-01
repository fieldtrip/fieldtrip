function state = dss_check_adaptivity(state)
% Check adaptivity of alpha, beta & gamma

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

state.adapt_alpha=0;
if isfield(state, 'alphaf')
    p = feval(state.alphaf.h);
    if p.adaptive==1
        state.adapt_alpha=1;
        dss_message(state, 1, 'Adaptive alpha\n');
    else
        [state.alphaf.params, state.alpha] = feval(state.alphaf.h, state.alphaf.params, state);
        dss_message(state, 1, sprintf('Constant alpha: %d\n', state.alpha));
    end
else dss_message(state, 1, sprintf('Fixed alpha: %d\n', state.alpha));
end

state.adapt_beta=0;
if isfield(state, 'betaf')
    p = feval(state.betaf.h);
    if p.adaptive==1
        state.adapt_beta=1;
        dss_message(state, 1, 'Adaptive beta\n');
    else
        [state.betaf.params, state.beta] = feval(state.betaf.h, state.betaf.params, state);
        dss_message(state, 1, sprintf('Constant beta: %d\n', state.beta));
    end
else dss_message(state, 1, sprintf('Fixed beta: %d\n', state.beta));
end

state.adapt_gamma=0;
if isfield(state, 'gammaf')
    p = feval(state.gammaf.h);
    if p.adaptive==1
        state.adapt_gamma=1;
        dss_message(state, 1, 'Adaptive gamma\n');
    else
        [state.gammaf.params, state.gamma] = feval(state.gammaf.h, state.gammaf.params, state);
        dss_message(state, 1, sprintf('Constant gamma: %d\n', state.gamma));
    end
else dss_message(state, 1, sprintf('Fixed gamma: %d\n', state.gamma));
end


% -- Initial values for alpha, beta & gamma
%for f = {'alpha', 'beta', 'gamma'}
%  fname = [f{1} 'f'];
%  if isfield(state, 'alphaf')
%    dss_message(state, 1, [fname ' is @' func2str(state.(fname).h) '.\n']);
%    [state.(f{1}), state.(fname).params] = ...
%	feval(state.(fname).h, state.(fname).params, state);
%  else
%    dss_message(state, 1, [fname ' is not defined.\n']);
%  end
%end
