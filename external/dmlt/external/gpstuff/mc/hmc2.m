function [samples, energies, diagn] = hmc2(f, x, opt, gradf, varargin)
%HMC2   Hybrid Monte Carlo sampling.
%
%       Description
%       SAMPLES = HMC2(F, X, OPTIONS, GRADF) uses a  hybrid Monte Carlo
%       algorithm to sample from the distribution P ~ EXP(-F), where F is the
%       first argument to HMC2. The Markov chain starts at the point X, and
%       the function GRADF is the gradient of the `energy' function F.
%
%       HMC2(F, X, OPTIONS, GRADF, P1, P2, ...) allows additional arguments to
%       be passed to F() and GRADF().
%
%       [SAMPLES, ENERGIES, DIAGN] = HMC2(F, X, OPTIONS, GRADF) also returns a
%       log of the energy values (i.e. negative log probabilities) for the
%       samples in ENERGIES and DIAGN, a structure containing diagnostic
%       information (position, momentum and acceptance threshold) for each
%       step of the chain in DIAGN.POS, DIAGN.MOM and DIAGN.ACC respectively.
%       All candidate states (including rejected ones) are stored in
%       DIAGN.POS. The DIAGN structure contains fields: 
%
%       pos
%        the position vectors of the dynamic process
%       mom
%        the momentum vectors of the dynamic process
%       acc
%        the acceptance thresholds
%       rej
%        the number of rejections
%       stp
%        the step size vectors
%
%       S = HMC2('STATE') returns a state structure that contains
%       the used random stream and its state and the momentum of
%       the dynamic process. These are contained in fields stream,
%       streamstate and mom respectively. The momentum state is
%       only used for a persistent momentum update.
%
%       HMC2('STATE', S) resets the state to S. If S is an integer,
%       then it is used as a seed for the random stream and the
%       momentum variable is randomised. If S is a structure
%       returned by HMC2('STATE') then it resets the random stream
%       to exactly the same state.
%
%       See HMC2_OPT for the optional parameters in the OPTIONS structure.
%
%       See also
%       HMC2_OPT, METROP
%

%       Copyright (c) 1996-1998 Christopher M Bishop, Ian T Nabney
%       Copyright (c) 1998-2000,2010 Aki Vehtari

%  The portion of the HMC algorithm involving "windows" is derived
%  from the C code for this function included in the Software for
%  Flexible Bayesian Modeling written by Radford Neal
%  <http://www.cs.toronto.edu/~radford/fbm.software.html>. 

% Global variable to store state of momentum variables: set by set_state
% Used to initialize variable if set

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

global HMC_MOM

%Return state or set state to given state.
if nargin <= 2
  if ~strcmp(f, 'state')
    error('Unknown argument to hmc2');
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

% Set empty options to default values
opt=hmc2_opt(opt);
% Reference to structures is much slower, so...
opt_nsamples=opt.nsamples;
opt_window=opt.window;
opt_steps=opt.steps;
opt_display = opt.display;
opt_persistence = opt.persistence;

if opt_persistence
  alpha = opt.decay;
  salpha = sqrt(1-alpha*alpha);
end

% Stepsizes, varargin gives the opt.stepsf arguments (net, x ,y)
% where x is input data and y is a target data.
if ~isempty(opt.stepsf)
  epsilon = feval(opt.stepsf,varargin{:}).*opt.stepadj;
else
  epsilon = opt.stepadj;
end

% Force x to be a row vector
x = x(:)';
nparams = length(x);

%Set up strings for evaluating potential function and its gradient.
f = fcnchk(f, length(varargin));
gradf = fcnchk(gradf, length(varargin));

% Check the gradient evaluation.
if (opt.checkgrad)
  % Check gradients
  feval('gradcheck', x, f, gradf, varargin{:});
end

% Matrix of returned samples.
samples = zeros(opt_nsamples, nparams);
% Return energies?
if nargout >= 2
  en_save = 1;
  energies = zeros(opt_nsamples, 1);
else
  en_save = 0;
end
% Return diagnostics?
if nargout >= 3
  diagnostics = 1;
  diagn_pos = zeros(opt_nsamples, nparams);
  diagn_mom = zeros(opt_nsamples, nparams);
  diagn_acc = zeros(opt_nsamples, 1);
else
  diagnostics = 0;
end

if (~opt_persistence | isempty(HMC_MOM) | nparams ~= length(HMC_MOM))
  % Initialise momenta at random
  p = randn(1, nparams);
else
  % Initialise momenta from stored state
  p = HMC_MOM;
end

% Evaluate starting energy.
E = feval(f, x, varargin{:});

% Main loop.
nreject = 0;         %number of rejected samples
window_offset=0;     %window offset initialised to zero
k = - opt.nomit + 1;       %nomit samples are omitted, so we store 
while k <= opt_nsamples    %samples from k>0

  % Store starting position and momenta
  xold = x;
  pold = p;
  % Recalculate Hamiltonian as momenta have changed
  Eold = E;
  Hold = E + 0.5*(p*p');

  % Decide on window offset, if windowed HMC is used
  if opt_window>1
    window_offset=fix(opt_window*rand(1));
  end

  have_rej = 0;
  have_acc = 0;
  n = window_offset;
  dir = -1;           %the default value for dir (=direction)
                      %assumes that windowing is used
  while (dir==-1 | n~=opt_steps)

    %if windowing is not used or we have allready taken
    %window_offset steps backwards...
    if (dir==-1 & n==0)
      % Restore, next state should be original start state.
      if window_offset > 0
	x = xold;
	p = pold;
	n = window_offset;
      end
      %set dir for forward steps
      E = Eold;
      H = Hold;
      dir = 1;
      stps = dir;
    else
      if (n*dir+1<opt_window | n>(opt_steps-opt_window))
	% State in the accept and/or reject window.
	stps = dir;
      else
	% State not in the accept and/or reject window. 
	stps = opt_steps-2*(opt_window-1);
      end
      % First half-step of leapfrog.
      p = p - dir*0.5*epsilon.*feval(gradf, x, varargin{:});
      x = x + dir*epsilon.*p;
      % Full leapfrog steps.
      for m = 1:(abs(stps)-1)
	p = p - dir*epsilon.*feval(gradf, x, varargin{:});
	x = x + dir*epsilon.*p;
      end
      % Final half-step of leapfrog.
      p = p - dir*0.5*epsilon.*feval(gradf, x, varargin{:});

      E = feval(f, x, varargin{:});
      H = E + 0.5*(p*p');
      
      n=n+stps;
    end


    if (opt_window~=opt_steps+1 & n<opt_window)
      % Account for state in reject window.  Reject window can be
      % ignored if windows consist of the entire trajectory.
      if ~have_rej
	rej_free_energy = H;
      else
	rej_free_energy = -addlogs(-rej_free_energy, -H);
      end
      if (~have_rej | rand(1) < exp(rej_free_energy-H));
	E_rej=E;
	x_rej=x;
	p_rej=p;
        have_rej = 1;
      end
    end

    if (n>(opt_steps-opt_window))
      % Account for state in the accept window.
      if ~have_acc
	acc_free_energy = H;
      else
	acc_free_energy = -addlogs(-acc_free_energy, -H);
      end
      if (~have_acc | rand(1) < exp(acc_free_energy-H))
	E_acc=E;
	x_acc=x;
	p_acc=p;
	have_acc = 1;
      end
    end
  end

  % Acceptance threshold.
  a = exp(rej_free_energy - acc_free_energy);

  if (diagnostics & k > 0)
    diagn_pos(k,:) = x_acc;
    diagn_mom(k,:) = p_acc;
    diagn_acc(k,:) = a;
  end
  if (opt_display > 1)
    fprintf(1, 'New position is\n');
    disp(x);
  end

  % Take new state from the appropriate window.
  if a > rand(1)
    % Accept 
    E=E_acc;
    x=x_acc;
    p=-p_acc; % Reverse momenta
    if (opt_display > 0)
      fprintf(1, 'Finished step %4d  Threshold: %g\n', k, a);
    end
  else
    % Reject
    if k > 0
      nreject = nreject + 1;
    end
    E=E_rej;
    x=x_rej;
    p=p_rej;
    if (opt_display > 0)
      fprintf(1, '  Sample rejected %4d.  Threshold: %g\n', k, a);
    end
  end

  if k > 0
    % Store sample
    samples(k,:) = x;
    if en_save
      % Store energy
      energies(k) = E;
    end
  end

  % Set momenta for next iteration
  if opt_persistence
    % Reverse momenta
    p = -p;
    % Adjust momenta by a small random amount
    p = alpha.*p + salpha.*randn(1, nparams);
  else
    % Replace all momenta
    p = randn(1, nparams);
  end

  k = k + 1;
end

if (opt_display > 0)
  fprintf(1, '\nFraction of samples rejected:  %g\n', ...
    nreject/(opt_nsamples));
end

% Store diagnostics
if diagnostics
  diagn.pos = diagn_pos;   %positions matrix
  diagn.mom = diagn_mom;   %momentum matrix
  diagn.acc = diagn_acc;   %acceptance treshold matrix
  diagn.rej = nreject/(opt_nsamples);   %rejection rate
  diagn.stps = epsilon;    %stepsize vector
end

% Store final momentum value in global so that it can be retrieved later
if opt_persistence
  HMC_MOM = p;
else
  HMC_MOM = [];
end
return


function state = get_state(f)
%GET_STATE           Return complete state of sampler 
%                    (including momentum)

global HMC_MOM
state.stream=setrandstream();
state.streamstate = state.stream.State;
state.mom = HMC_MOM;
return


function set_state(f, x)
%SET_STATE    Set complete state of sample
%
%          Description
%          Set complete state of sampler (including momentum) 
%          or just set randn and rand with integer argument.
global HMC_MOM
if isnumeric(x)
  setrandstream(x);
else
  if ~isstruct(x)
    error('Second argument to hmc must be number or state structure');
  end
  if (~isfield(x, 'stream') | ~isfield(x, 'streamstate') ...
      | ~isfield(x, 'mom'))
    error('Second argument to hmc must contain correct fields')
  end
  setrandstream(x.stream);
  x.State=x.streamstate;
  HMC_MOM = x.mom;
end
return


function c=addlogs(a,b)
%ADDLOGS(A,B)  Add numbers represented by their logarithms.
%
%        Description
%        Add numbers represented by their logarithms.
%        Computes log(exp(a)+exp(b)) in such a fashion that it 
%        works even when a and b have large magnitude.
if a>b
  c = a + log(1+exp(b-a));
else
  c = b + log(1+exp(a-b));
end
