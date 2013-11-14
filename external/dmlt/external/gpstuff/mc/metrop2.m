function [samples, energies, diagn] = metrop2(f, x, opt, gradf, varargin)
%METROP2	Markov Chain Monte Carlo sampling with Metropolis algorithm.
%
%	Description
%	SAMPLES = METROP(F, X, OPT) uses the Metropolis algorithm to
%	sample from the distribution P ~ EXP(-F), where F is the first
%	argument to METROP.   The Markov chain starts at the point X and each
%	candidate state is picked from a Gaussian proposal distribution and
%	accepted or rejected according to the Metropolis criterion.
%
%	SAMPLES = METROP(F, X, OPT, [], P1, P2, ...) allows additional
%	arguments to be passed to F().  The fourth argument is ignored, but
%	is included for compatibility with HMC and the optimizers.
%
%	[SAMPLES, ENERGIES, DIAGN] = METROP(F, X, OPT) also returns a log
%	of the energy values (i.e. negative log probabilities) for the
%	samples in ENERGIES and DIAGN, a structure containing diagnostic
%	information (position and acceptance threshold) for each step of the
%	chain in DIAGN.POS and DIAGN.ACC respectively.  All candidate states
%	(including rejected ones) are stored in DIAGN.POS.
%
%	S = METROP('STATE') returns a state structure that contains the state
%	of the two random number generators RAND and RANDN. These are
%	contained in fields randstate,  randnstate.
%
%	METROP('STATE', S) resets the state to S.  If S is an integer, then
%	it is passed to RAND and RANDN. If S is a structure returned by
%	METROP('STATE') then it resets the generator to exactly the same
%	state.
%
%       See METROP2_OPT for the optional parameters in the OPTIONS
%       structure.
%
%	See also
%	HMC, METROP2_OPT
%

%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)
%	Copyright (c) 1998-2000 Aki Vehtari 

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Global variable to store state of momentum variables: set by set_state
% Used to initialize variable if set
global HMC_MOM
if nargin <= 2
  if ~strcmp(f, 'state')
    error('Unknown argument to metrop2');
  end
  switch nargin
   case 1
    samples = get_state(f);
    return;
   case 2
    set_state(f, x);
    return;
  end
end

% Set empty omptions to default values
opt=metrop2_opt(opt);
% Refrence to structures is much slower, so...
opt_nsamples=opt.nsamples;
opt_display =opt.display;
stddev = opt.stddev;


% Set up string for evaluating potential function.
%f = fcnchk(f, length(varargin));

nparams = length(x);
samples = zeros(opt_nsamples, nparams);		% Matrix of returned samples.
if nargout >= 2
  en_save = 1;
  energies = zeros(opt_nsamples, 1);
else
  en_save = 0;
end
if nargout >= 3
  diagnostics = 1;
  diagn_pos = zeros(opt_nsamples, nparams);
  diagn_acc = zeros(opt_nsamples, 1);
else
  diagnostics = 0;
end

% Main loop.
k = - opt.nomit + 1;
Eold = f(x, varargin{:});	% Evaluate starting energy.
nreject = 0;				% Initialise count of rejected states.
while k <= opt_nsamples

  xold = x;
  % Sample a new point from the proposal distribution
  x = xold + randn(1, nparams)*stddev;

  % Now apply Metropolis algorithm.
  Enew = f(x, varargin{:});	% Evaluate new energy.
  a = exp(Eold - Enew);			% Acceptance threshold.
  if (diagnostics & k > 0)
    diagn_pos(k,:) = x;
    diagn_acc(k,:) = a;
  end
  if (opt_display > 1)
    fprintf(1, 'New position is\n');
    disp(x);
  end

  if a > rand(1)	% Accept the new state.
    Eold = Enew;
    if (opt_display > 0)
      fprintf(1, 'Finished step %4d  Threshold: %g\n', k, a);
    end
  else			% Reject the new state
    if k > 0
      nreject = nreject + 1;
    end
    x = xold;	% Reset position 
    if (opt_display > 0)
      fprintf(1, '  Sample rejected %4d.  Threshold: %g\n', k, a);
    end
  end
  if k > 0
    samples(k,:) = x;			% Store sample.
    if en_save 
      energies(k) = Eold;		% Store energy.
    end
  end
  k = k + 1;
end

if (opt_display > 0)
  fprintf(1, '\nFraction of samples rejected:  %g\n', ...
          nreject/(opt_nsamples));
end

if diagnostics
  diagn.pos = diagn_pos;
  diagn.acc = diagn_acc;
end

% Return complete state of the sampler.
function state = get_state(f)

state.randstate = rand('state');
state.randnstate = randn('state');
return

% Set state of sampler, either from full state, or with an integer
function set_state(f, x)

if isnumeric(x)
  rand('state', x);
  randn('state', x);
else
  if ~isstruct(x)
    error('Second argument to metrop must be number or state structure');
  end
  if (~isfield(x, 'randstate') | ~isfield(x, 'randnstate'))
    error('Second argument to metrop must contain correct fields')
  end
  rand('state', x.randstate);
  randn('state', x.randnstate);
end
return
