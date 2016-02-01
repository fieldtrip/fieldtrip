function [samples,energies,diagn] = sls(f, x, opt, gradf, varargin)
%SLS  Markov Chain Monte Carlo sampling using Slice Sampling
%
%  Description
%    SAMPLES = SLS(F, X, OPTIONS) uses slice sampling to sample
%      from the distribution P ~ EXP(-F), where F is the first
%      argument to SLS. Markov chain starts from point X and the
%      sampling from multivariate distribution is implemented by
%      sampling each variable at a time either using overrelaxation
%      or not. See SLS_OPT for details. A simple multivariate scheme
%      using hyperrectangles is utilized when method is defined 'multi'.
%
%    SAMPLES = SLS(F, X, OPTIONS, [], P1, P2, ...) allows additional
%      arguments to be passed to F(). The fourth argument is ignored,
%      but included for compatibility with HMC and the optimisers.
%
%    [SAMPLES, ENERGIES, DIAGN] = SLS(F, X, OPTIONS) Returns some additional
%      diagnostics for the values in SAMPLES and ENERGIES.
%
%  See SLS_OPT for the optional parameters in the OPTIONS structure.
%
%  See also
%    METROP2, HMC2, SLS_OPT

%  Based on "Slice Sampling" by Radford M. Neal in "The Annals of Statistics"
%  2003, Vol. 31, No. 3, 705-767, (c) Institute of Mathematical Statistics, 2003
%  Thompson & Neal (2010) Covariance-Adaptive Slice Sampling. Technical
%  Report No. 1002, Department of Statistic, University of Toronto
%
%  Copyright (c) Toni Auranen, 2003-2006
%  Copyright (c) Ville Tolvanen, 2012

% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

%  Version 1.01, 21/1/2004, TA
%  Version 1.02, 28/1/2004, TA
%   - fixed the limit-checks for stepping_out and doubling
%  Version 1.02b, 26/2/2004, TA
%   - changed some display-settings
%  Version 1.03, 9/3/2004, TA
%   - changed some display-settings
%   - changed the variable initialization
%   - optimized the number of fevals
%  Version 1.04, 24/3/2004, TA
%   - overrelaxation separately for each variable in multivariate case
%   - added maxiter parameter for shrinkage
%   - added some argument checks
%  Version 1.04b, 29/3/2004, TA
%   - minor fixes
%  Version 1.05, 7/4/2004, TA
%   - minor bug fixes
%   - added nomit-option
%  Version 1.06, 22/4/2004, TA
%   - added unimodality shortcut for stepping-out and doubling
%   - optimized the number of fevals in doubling
%  Version 1.06b, 27/4/2004, TA
%   - fixed some bugs
%  Version 1.7, 15/3/2005, TA
%   - added the hyperrectangle multivariate sampling


% Start timing and construct a function handle from the function name string
% (Timing is off, and function handles are left for the user)
%t = cputime;
%f = str2func(f);

% Set empty options to default values
opt = sls_opt(opt);
%if opt.display, disp(opt); end
if opt.display == 1
  opt.display = 2; % verbose
elseif opt.display == 2
  opt.display = 1; % all
end

% Forces x to be a row vector
x = x(:)';

% Set up some variables
nparams = length(x);
samples = zeros(opt.nsamples,nparams);
if nargout >= 2
  save_energies = 1;
  energies = zeros(opt.nsamples,1);
else
  save_energies = 0;
end
if nargout >= 3
  save_diagnostics = 1;
else
  save_diagnostics = 0;
end
if nparams == 1
  multivariate = 0;
  if strcmp(opt.method,'multi')
    opt.method = 'stepping';
  end
end
if nparams > 1
  multivariate = 1;
end
rej = 0;
rej_step = 0;
rej_old = 0;
x_0 = x;

umodal = opt.unimodal;
nomit = opt.nomit;
nsamples = opt.nsamples;
display_info = opt.display;
method = opt.method;
overrelaxation = opt.overrelaxation;
overrelaxation_info = ~isempty(find(overrelaxation));
w = opt.wsize;
maxiter = opt.maxiter;
m = opt.mlimit;
p = opt.plimit;
a = opt.alimit;
mmin = opt.mmlimits(1,:);
mmax = opt.mmlimits(2,:);

if multivariate
  if length(w) == 1
    w = w.*ones(1,nparams);
  end
  if length(m) == 1
    m = m.*ones(1,nparams);
  end
  if length(p) == 1
    p = p.*ones(1,nparams);
  end
  if length(overrelaxation) == 1
    overrelaxation = overrelaxation.*ones(1,nparams);
  end
  if length(a) == 1
    a = a.*ones(1,nparams);
  end
  if length(mmin) == 1
    mmin = mmin.*ones(1,nparams);
  end
  if length(mmax) == 1
    mmax = mmax.*ones(1,nparams);
  end  
end
if overrelaxation_info
  nparams_or = length(find(overrelaxation));
end
if ~isempty(find(w<=0))
  error('Parameter ''wsize'' must be positive.');
end
if (strcmp(method,'stepping') || strcmp(method,'doubling')) && isempty(find(mmax-mmin>2*w))
  error('Check parameter ''mmlimits''. The interval is too small in comparison to parameter ''wsize''.');
end
if strcmp(method,'stepping') && ~isempty(find(m<1))
  error('Parameter ''mlimit'' must be >0.');
end
if overrelaxation_info && ~isempty(find(a<1))
  error('Parameter ''alimit'' must be >0.');
end
if strcmp(method,'doubling') && ~isempty(find(p<1))
  error('Parameter ''plimit'' must be >0.');
end

ind_umodal = 0;
j = 0;
y_new = -f(x_0,varargin{:});

% The main loop of slice sampling
for i = 1-nomit:1:nsamples
  switch method
   % Slice covariance matching from Thompson & Neal (2010)
   case 'covmatch' 
    theta = 1;
    np = length(x_0);
    M = y_new;
    ee = exprnd(1);
    ytilde0 = M-ee;
    x_0 = x_0';
    sigma = opt.sigma;
    R = 1./sigma*eye(np);
    F = R;
    cbarstar = 0;
    
    while 1
      z = mvnrnd(zeros(1,np), eye(np))';
      c = x_0 + F\z;
      cbarstar = cbarstar + F'*(F*c);
      cbar = R\(R'\cbarstar);
      z = mvnrnd(zeros(1,np), eye(np))';
      x_prop = cbar + R\z;
      y_new = -f(x_prop',varargin{:});
      if y_new > ytilde0
        % Accept proposal
        break;
      end
      G = -gradf(x_prop', varargin{:})';
      gr = G/norm(G);
      delta = norm(x_prop - c);
      u = x_prop + delta*gr;
      lu = -f(u', varargin{:});
      kappa = -2/delta^2*(lu-y_new-delta*norm(G));
      lxu = 0.5*norm(G)^2/kappa + y_new;
      M = max(M, lxu);
      sigma2 = 2/3*(M-ytilde0)/kappa;
      alpha = max(0, 1/sigma2 - (1+theta)*gr'*(R'*(R*gr)));
      F = chol(theta*(R'*R) + alpha*(gr*gr'));
      R = chol((1+theta)*(R'*R) + alpha*(gr*gr'));
      
    end
    % Save sampling step and set up the new 'old' sample
    x_0 = x_prop;
    if i > 0
      samples(i,:) = x_prop;
    end
    
    % Save energies
    if save_energies && i > 0
      energies(i) = -y_new;
    end
    
   % Shrinking-Rank method from Thompson & Neal (2010)
   case 'shrnk'
    ytr = -log(rand) - y_new;
    k = 0;
    sigma(1) = opt.sigma;
    J = [];
    np = length(x_0);
    
    while 1
      k = k+1;
      c(k,:) = P(J, mvnrnd(x_0, sigma(k).^2*eye(np)));
      sigma2 = 1./(sum(1./sigma.^2));
      mu = sigma2*(sum(bsxfun(@times, 1./sigma'.^2, bsxfun(@minus, c,x_0)),1));
      x_prop = x_0 + P(J, mvnrnd(mu, sqrt(sigma2)*eye(np)));
      y_new = f(x_prop, varargin{:});
      if y_new < ytr
        % Accept proposal
        break;
      end
      gradient = gradf(x_prop, varargin{:});
      gstar = P(J, gradient);
      if size(J,2) < size(x,2)-1 && gstar*gradient'/(norm(gstar)*norm(gradient)) > cos(pi/3)
        J = [J gstar'/norm(gstar)];
        sigma(k+1) = sigma(k);
      else
        sigma(k+1) = 0.95*sigma(k);
      end
      
    end
    % Save sampling step and set up the new 'old' sample
    x_0 = x_prop;
    if i > 0
      samples(i,:) = x_prop;
    end
    
    % Save energies
    if save_energies && i > 0
      energies(i) = y_new;
    end
    
   % Multivariate rectangle sampling step
   case 'multi'
    x_new = x_0;
    y = y_new + log(rand(1));
    if isinf(y)
      x_new = mmin + (mmax-mmin).*rand(1,length(x_new));
      y_new = -f(x_new,varargin{:});
    else
      L = max(x_0 - w.*rand(1,length(x_0)),mmin);
      R = min(L + w,mmax);
      x_new = L + rand(1,length(x_new)).*(R-L);
      y_new = -f(x_new,varargin{:});  
      while y >= y_new
        L(x_new < x_0) = x_new(x_new < x_0);
        R(x_new >= x_0) = x_new(x_new >= x_0);
        x_new = L + rand(1,length(x_new)).*(R-L);
        y_new = -f(x_new,varargin{:});        
      end % while
    end % isinf(y)
    
    % Save sampling step and set up the new 'old' sample
    x_0 = x_new;
    if i > 0
      samples(i,:) = x_new;
    end
    
    % Save energies
    if save_energies && i > 0
      energies(i) = -y_new;
    end
    
    % Display energy information
    if display_info == 1
      fprintf('Finished multi-step %4d  Energy: %g\n',i,-y_new);      
    end
        
   case 'multimm'
    x_new = x_0;
    y = y_new + log(rand(1));
    if isinf(y)
      x_new = mmin + (mmax-mmin).*rand(1,length(x_new));
      y_new = -f(x_new,varargin{:});
    else
      L = mmin;
      R = mmax;
      x_new = L + rand(1,length(x_new)).*(R-L);
      y_new = -f(x_new,varargin{:});  
      while y >= y_new
        L(x_new < x_0) = x_new(x_new < x_0);
        R(x_new >= x_0) = x_new(x_new >= x_0);
        x_new = L + rand(1,length(x_new)).*(R-L);
        y_new = -f(x_new,varargin{:});        
      end % while
    end % isinf(y)
    
    % Save sampling step and set up the new 'old' sample
    x_0 = x_new;
    if i > 0
      samples(i,:) = x_new;
    end
    
    % Save energies
    if save_energies && i > 0
      energies(i) = -y_new;
    end
    
    % Display energy information
    if display_info == 1
      fprintf('Finished multimm-step %4d  Energy: %g\n',i,-y_new);      
    end
        
   % Other sampling steps
   otherwise
    ind_umodal = ind_umodal + 1;
    x_new = x_0;
    for j = 1:nparams
      y = y_new + log(rand(1));
      if isinf(y)
        x_new(j) = mmin(j) + (mmax(j)-mmin(j)).*rand;
        y_new = -f(x_new,varargin{:});
      else
        L = x_new;
        R = x_new;
        switch method
         case 'stepping'
          [L, R] = stepping_out(f,y,x_new,L,R,w,m,j,mmin,mmax,display_info,umodal,varargin{:});
         case 'doubling'
          [L, R] = doubling(f,y,x_new,L,R,w,p,j,mmin,mmax,display_info,umodal,varargin{:});
         case 'minmax'
          L(j) = mmin(j);
          R(j) = mmax(j);
        end % switch
        if overrelaxation(j)
          [x_new, y_new, rej_step, rej_old] = bisection(f,y,x_new,L,R,w,a,rej_step,j,umodal,varargin{:});        
        else
          [x_new, y_new] = shrinkage(f,y,x_new,w,L,R,method,j,maxiter,umodal,varargin{:});
        end % if overrelaxation
        if umodal % adjust the slice if the distribution is known to be unimodal
          w(j) = (w(j)*ind_umodal + abs(x_0(j)-x_new(j)))/(ind_umodal+1);
        end % if umodal
      end % if isinf(y)
    end % j:nparams
    if overrelaxation_info & multivariate
      rej = rej + rej_step/nparams_or;
    elseif overrelaxation_info & ~multivariate
      rej = rej + rej_step;
    end
    
    % Save sampling step and set up the new 'old' sample
    x_0 = x_new;
    if i > 0
      samples(i,:) = x_new;
    end
    
    % Save energies
    if save_energies && i > 0
      energies(i) = -y_new;
    end
    
    % Display information and keep track of rejections (overrelaxation)
    if display_info == 1
      if ~multivariate && overrelaxation_info && rej_old
        fprintf('    Sample %4d rejected (overrelaxation).\n',i);
        rej_old = 0;
        rej_step = 0;
      elseif multivariate && overrelaxation_info
        fprintf('Finished step %4d (RR: %1.1f, %d/%d)  Energy: %g\n',i,100*rej_step/nparams_or,nparams_or,nparams,-y_new);
        rej_step = 0;
        rej_old = 0;
      else
        fprintf('Finished step %4d  Energy: %g\n',i,-y_new);      
      end
    else
      rej_old = 0;
      rej_step = 0;
    end
  
  end % switch
end % i:nsamples

% Save diagnostics
if save_diagnostics
  diagn.opt = opt;
end

% Display rejection information after slice sampling is complete (overrelaxation)
if overrelaxation_info && nparams == 1 && display_info == 1
  fprintf('\nRejected samples due to overrelaxation (percentage): %1.1f\n',100*rej/nsamples);
elseif overrelaxation_info && nparams > 1 && display_info == 1
  fprintf('\nAverage rejections per step due to overrelaxation (percentage): %1.1f\n',100*rej/nsamples);
end

% Display the elapsed time
%if display_info == 1
%  if (cputime-t)/60 < 4
%    fprintf('\nElapsed cputime (seconds): %1.1f\n\n',cputime-t);
%  else
%    fprintf('\nElapsed cputime (minutes): %1.1f\n\n',(cputime-t)/60);
%  end
%end

%disp(w);

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function [x_new, y_new, rej, rej_old] = bisection(f,y,x_0,L,R,w,a,rej,j,um,varargin);
%function [x_new, y_new, rej, rej_old] = bisection(f,y,x_0,L,R,w,a,rej,j,um,varargin);
%
% Bisection for overrelaxation (stepping-out needs to be used)

x_new = x_0;
M = (L + R) / 2;
l = L;
r = R;
q = w(j);
s = a(j);
if (R(j) - L(j)) < 1.1*w(j)
  while 1
    M(j) = (l(j) + r(j))/2;
    if s == 0 || y < -f(M,varargin{:})
      break;
    end
    if x_0(j) > M(j)
      l(j) = M(j);
    else
      r(j) = M(j);
    end
    s = s - 1;
    q = q / 2;
  end % while
end % if
ll = l;
rr = r;
while s > 0
  s = s - 1;
  q = q / 2;
  tmp_ll = ll;
  tmp_ll(j) = tmp_ll(j) + q;
  tmp_rr = rr;
  tmp_rr(j) = tmp_rr(j) - q;
  if y >= -f((tmp_ll),varargin{:})
    ll(j) = ll(j) + q;
  end
  if y >= -f((tmp_rr),varargin{:})
    rr(j) = rr(j) - q;
  end
end % while
x_new(j) = ll(j) + rr(j) - x_0(j);
y_new = -f(x_new,varargin{:});
if x_new(j) < l(j) || x_new(j) > r(j) || y >= y_new
  x_new(j) = x_0(j);
  rej = rej + 1;
  rej_old = 1;
  y_new = y;
else
  rej_old = 0;
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function [x_new, y_new] = shrinkage(f,y,x_0,w,L,R,method,j,maxiter,um,varargin);
%function [x_new, y_new] = shrinkage(f,y,x_0,w,L,R,method,j,maxiter,um,varargin);
%    
% Shrinkage with acceptance-check for doubling scheme
% - acceptance-check is skipped if the distribution is defined
%   to be unimodal by the user

iter = 0;
x_new = x_0;
l = L(j);
r = R(j);
while 1
  x_new(j) = l + (r-l).*rand;
  if strcmp(method,'doubling')
    y_new = -f(x_new,varargin{:});
    if y < y_new && (um || accept(f,y,x_0,x_new,w,L,R,j,varargin{:}))
      break;
    end
  else
    y_new = -f(x_new,varargin{:});
    if y < y_new
      break;
      break;
    end
  end % if strcmp
  if x_new(j) < x_0(j)
    l = x_new(j);
  else
    r = x_new(j);
  end % if
  iter = iter + 1;
  if iter > maxiter
    fprintf('Maximum number (%d) of iterations reached for parameter %d during shrinkage.\n',maxiter,j);
    if strcmp(method,'minmax')
      error('Check function F, decrease the interval ''mmlimits'' or increase the value of ''maxiter''.');
    else
      error('Check function F or increase the value of ''maxiter''.');
    end
  end
end % while

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function [L,R] = stepping_out(f,y,x_0,L,R,w,m,j,mmin,mmax,di,um,varargin);
%function [L,R] = stepping_out(f,y,x_0,L,R,w,m,j,mmin,mmax,di,um,varargin);
%
% Stepping-out procedure

if um % if the user defines the distribution to be unimodal
  L(j) = x_0(j) - w(j).*rand;
  if L(j) < mmin(j)
    L(j) = mmin(j);
    if di
      fprintf('Underflow! (L:%d)\n',j);
    end
  end
  R(j) = L(j) + w(j);
  if R(j) > mmax(j)
    R(j) = mmax(j);
    if di
      fprintf('Overflow! (R:%d)\n',j);
    end
  end
  while y < -f(L,varargin{:})
    L(j) = L(j) - w(j);
    if L(j) < mmin(j)
      L(j) = mmin(j);
      if di
        fprintf('Underflow! (L:%d)\n',j);
      end
      break;
    end
  end
  while y < -f(R,varargin{:})
    R(j) = R(j) + w(j);
    if R(j) > mmax(j)
      R(j) = mmax(j);
      if di
        fprintf('Overflow! (R:%d)\n',j);
      end
      break;
    end
  end
else % if the distribution is not defined to be unimodal
  L(j) = x_0(j) - w(j).*rand;
  J = floor(m(j).*rand);
  if L(j) < mmin(j)
    L(j) = mmin(j);
    if di
      fprintf('Underflow! (L:%d)\n',j);
    end
    J = 0;
  end
  R(j) = L(j) + w(j);
  K = (m(j)-1) - J;
  if R(j) > mmax(j)
    R(j) = mmax(j);
    if di
      fprintf('Overflow! (R:%d)\n',j);
    end
    K = 0;
  end
  while J > 0 && y < -f(L,varargin{:})
    L(j) = L(j) - w(j);
    if L(j) < mmin(j)
      L(j) = mmin(j);
      if di
        fprintf('Underflow! (L:%d)\n',j);
      end
      break;
    end
    J = J - 1;
  end
  while K > 0 && y < -f(R,varargin{:})
    R(j) = R(j) + w(j);
    if R(j) > mmax(j)
      R(j) = mmax(j);
      if di
        fprintf('Overflow! (R:%d)\n',j);
      end
      break;
    end
    K = K - 1;
  end
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function [L,R] = doubling(f,y,x_0,L,R,w,p,j,mmin,mmax,di,um,varargin);
%function [L,R] = doubling(f,y,x_0,L,R,w,p,j,mmin,mmax,di,um,varargin);
%
% Doubling scheme for slice sampling

if um % if the user defines the distribution to be unimodal
  L(j) = x_0(j) - w(j).*rand;
  if L(j) < mmin(j)
    L(j) = mmin(j);
    if di
      fprintf('Underflow! (L:%d)\n',j);
    end
    Ao = 1;
  else
    Ao = 0;
  end
  R(j) = L(j) + w(j);
  if R(j) > mmax(j)
    R(j) = mmax(j);
    if di
      fprintf('Overflow! (R:%d)\n',j);
    end
    Bo = 1;
  else
    Bo = 0;
  end
  AL = -f(L,varargin{:});
  AR = -f(R,varargin{:});
  while (Ao == 0 && y < AL) || (Bo == 0 && y < AR)
    if rand < 1/2
      L(j) = L(j) - (R(j)-L(j));
      if L(j) < mmin(j)
        L(j) = mmin(j);
        if di
          fprintf('Underflow! (L:%d)\n',j);
        end
        Ao = 1;
      else
        Ao = 0;
      end
      AL = -f(L,varargin{:});
    else
      R(j) = R(j) + (R(j)-L(j));
      if R(j) > mmax(j)
        R(j) = mmax(j);
        if di
          fprintf('Overflow! (R:%d)\n',j);
        end
        Bo = 1;
      else
        Bo = 0;
      end
      AR = -f(R,varargin{:});
    end
  end % while
else % if the distribution is not defined to be unimodal
  L(j) = x_0(j) - w(j).*rand;
  if L(j) < mmin(j)
    L(j) = mmin(j);
    if di
      fprintf('Underflow! (L:%d)\n',j);
    end
  end
  R(j) = L(j) + w(j);
  if R(j) > mmax(j)
    R(j) = mmax(j);
    if di
      fprintf('Overflow! (R:%d)\n',j);
    end
  end
  K = p(j);
  AL = -f(L,varargin{:});
  AR = -f(R,varargin{:});
  while K > 0 && (y < AL || y < AR)
    if rand < 1/2
      L(j) = L(j) - (R(j)-L(j));
      if L(j) < mmin(j)
        L(j) = mmin(j);
        if di
          fprintf('Underflow! (L:%d)\n',j);
        end
      end
      AL = -f(L,varargin{:});
    else
      R(j) = R(j) + (R(j)-L(j));
      if R(j) > mmax(j)
        R(j) = mmax(j);
        if di
          fprintf('Overflow! (R:%d)\n',j);
        end
      end
      AR = -f(R,varargin{:});
    end
    K = K - 1;
  end % while
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function out = accept(f,y,x_0,x_new,w,L,R,j,varargin)
%function out = accept(f,y,x_0,x_new,w,L,R,j,varargin)
%
% Acceptance check for doubling scheme

out = [];

l = L;
r = R;
d = 0;
while r(j)-l(j) > 1.1*w(j)
  m = (l(j)+r(j))/2;
  if (x_0(j) < m && x_new(j) >= m) || (x_0(j) >= m && x_new(j) < m)
    d = 1;
  end
  if x_new(j) < m
    r(j) = m;
  else
    l(j) = m;
  end
  if d && y >= -f(l,varargin{:}) && y >= -f(r,varargin{:})
    out = 0;
    break;
  end
end % while

if isempty(out)
  out = 1;
end;

function p = P(J,v)
if size(J,2) ~= 0
  p = v' - J*(J'*v');
  p = p';
else
  p = v;
end
