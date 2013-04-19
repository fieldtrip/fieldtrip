function [samples,samplesf,diagn] = surrogate_sls(f, x, opt, gp, xx, yy, z, varargin)
%SURROGATE_SLS  Markov Chain Monte Carlo sampling using Surrogate data Slice Sampling
%
%  Description
%    SAMPLES = SURROGATE_SLS(F, X, OPTIONS) uses slice sampling to sample
%      from the distribution P ~ EXP(-F), where F is the first
%      argument to SLS. Markov chain starts from point X and the
%      sampling from multivariate distribution is implemented by
%      sampling each variable at a time either using overrelaxation
%      or not. See SLS_OPT for details. A simple multivariate scheme
%      using hyperrectangles is utilized when method is defined 'multi'.
%
%    SAMPLES = SURROGATE_SLS(F, X, OPTIONS, [], P1, P2, ...) allows additional
%      arguments to be passed to F(). The fourth argument is ignored,
%      but included for compatibility with HMC and the optimisers.
%
%    [SAMPLES, ENERGIES, DIAGN] = SLS(F, X, OPTIONS) Returns some additional
%      diagnostics for the values in SAMPLES and ENERGIES.
%
%  See SSLS_OPT and SLS_OPT for the optional parameters in the OPTIONS structure.
%
%  See also
%    METROP2, HMC2, SLS_OPT, SLS

%  Based on "Slice Sampling" by Radford M. Neal in "The Annals of Statistics"
%  2003, Vol. 31, No. 3, 705-767, (c) Institute of Mathematical Statistics, 2003
%  "Slice sampling covariance hyperparameters of latent Gaussian models"
%  by Iain Murray and Ryan P. Adams, 2010, Arxiv preprint arXiv:1006.0868

%  Copyright (c) Toni Auranen, 2003-2006
%  Copyright (c) Ville Tolvanen, 2012


% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.


% Set empty options to default values
opt = ssls_opt(opt);
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
n = size(gp.latentValues,1);
samples = zeros(opt.nsamples,nparams);
samplesf = zeros(n,opt.nsamples);
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
f_0 = f;

ncf = length(gp.cf);
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
if isfield(opt, 'scale')
  scale = opt.scale;
else
  scale = 5;
end

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
% y_new = -f(x_0,varargin{:});
% S = diag(ones(size(f)));
% The main loop of slice sampling
for i = 1-nomit:1:nsamples
  
  % Treshold
  [y, tmp, eta, g] = getY(gp, xx, yy, z, f_0, x_0, []);
  
%   fprintf('\n')
%   fprintf('Treshold: %g\n',y);

  switch method
   % Multivariate rectangle sampling step
   case 'multi'
    x_new = x_0;
    
    L = max(x_0 - w.*rand(1,length(x_0)),mmin);
    R = min(L + w,mmax);
    x_new = L + rand(1,length(x_new)).*(R-L);
    
    [y_new, f_new] = getY(gp,xx,yy,z,[], x_new, eta, g);
    
    while y >= y_new
%       disp(y_new)
      L(x_new < x_0) = x_new(x_new < x_0);
      R(x_new >= x_0) = x_new(x_new >= x_0);
      if sum(abs(L-R))<1e-8
        error('BUG DETECTED: Shrunk to minimum position and still not acceptable.');
      end
      x_new = L + rand(1,length(x_new)).*(R-L);
      
      [y_new, f_new] = getY(gp, xx, yy, z, [], x_new, eta, g);

    end % while
    
    
    % Save sampling step and set up the new 'old' sample
    x_0 = x_new;
    f_0 = f_new;
    if i > 0
      samples(i,:) = x_new;
      latent_opt = esls(opt.latent_opt);
      gp = gp_unpak(gp, x_new);
      for ii=1:opt.fsamples-1
        f_new = esls(f_new, latent_opt, gp, xx, yy, z);
      end
      samplesf(:,i) = f_new;
      f_0 = samplesf(:,end);
    end
    
    % Save energies
    %     if save_energies && i > 0
    %       energies(i) = -y_new;
    %     end
    
    % Display energy information
    if display_info == 1
      fprintf('Finished multi-step %4d  Energy: %g\n',i,-y_new);
    end
    
   case 'multimm'
     x_new = x_0;
     %     if isinf(y)
     %       x_new = mmin + (mmax-mmin).*rand(1,length(x_new));
     %       y_new = -f(x_new,varargin{:});
     %     else
     L = mmin;
     R = mmax;
     x_new = L + rand(1,length(x_new)).*(R-L);
     [y_new, f_new] = getY(gp, xx, yy, z, [], x_new, eta, g);
     while y >= y_new
       L(x_new < x_0) = x_new(x_new < x_0);
       R(x_new >= x_0) = x_new(x_new >= x_0);
       x_new = L + rand(1,length(x_new)).*(R-L);
       [y_new, f_new] = getY(gp, xx, yy, z, [], x_new, eta, g);
     end % while
     %     end % isinf(y)
    
    % Save sampling step and set up the new 'old' sample
%     fprintf('Accepted: %g\n',y_new);
    x_0 = x_new;
    f_0 = f_new;
    if i > 0
      samples(i,:) = x_new;
      latent_opt = esls(opt.latent_opt);
      gp = gp_unpak(gp, x_new);
      for ii=1:opt.fsamples-1
        f_new = esls(f_new, latent_opt, gp, xx, yy, z);
      end
      samplesf(:,i) = f_new;
      f_0 = samplesf(:,end);
    end
    
    % Save energies
%     if save_energies && i > 0
%       energies(i) = -y_new;
%     end
%     
    % Display energy information
    if display_info == 1
      fprintf('Finished multimm-step %4d  Energy: %g\n',i,-y_new);      
    end
        
   % Other sampling steps
   otherwise
    ind_umodal = ind_umodal + 1;
    x_new = x_0;
    f_new = f_0;
    for j = 1:nparams
      L = x_new;
      R = x_new;
      switch method
        case 'stepping'
          [L, R] = stepping_out(f_new,y,x_new,L,R,w,m,j,mmin,mmax,display_info,umodal,xx,yy,gp,z,eta,g,varargin{:});
        case 'doubling'
          [L, R] = doubling(f_new,y,x_new,L,R,w,m,j,mmin,mmax,display_info,umodal,xx,yy,gp,z,eta,g,varargin{:});
        case 'minmax'
          L(j) = mmin(j);
          R(j) = mmax(j);
        otherwise
          error('unknown method');
      end % switch
      if overrelaxation(j)
        [x_new, f_new, rej_step, rej_old, y_new] = bisection(f_new,y,x_new,L,R,w,a,rej_step,j,umodal,xx,yy,gp,z,eta,g);
      else
        [x_new, f_new] = shrinkage(f_new,y,x_new,w,L,R,method,j,maxiter,umodal,xx,yy,gp,z,eta,g);
      end % if overrelaxation
      if umodal % adjust the slice if the distribution is known to be unimodal
        w(j) = (w(j)*ind_umodal + abs(x_0(j)-x_new(j)))/(ind_umodal+1);
      end % if umodal
      %       end % if isinf(y)
    end % j:nparams
    if overrelaxation_info & multivariate
      rej = rej + rej_step/nparams_or;
    elseif overrelaxation_info & ~multivariate
      rej = rej + rej_step;
    end
    
    % Save sampling step and set up the new 'old' sample
    x_0 = x_new;
    f_0 = f_new;
    if i > 0
      samples(i,:) = x_new;
      latent_opt = esls(opt.latent_opt);
      gp = gp_unpak(gp, x_new);
      for ii=1:opt.fsamples-1
        f_new = esls(f_new, latent_opt, gp, xx, yy, z);
      end
      samplesf(:,i) = f_new;
      f_0 = f_new;
    end
    
    % Save energies
%     if save_energies && i > 0
%       energies(i) = -y_new;
%     end
    
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

function [x_new, f_new, rej, rej_old, y_new] = bisection(f,y,x_0,L,R,w,a,rej,j,um,xx,yy,gp,z,eta,g);
%function [x_new, y_new, rej, rej_old] = bisection(f,y,x_0,L,R,w,a,rej,j,um,varargin);
%
% Bisection for overrelaxation (stepping-out needs to be used)

x_new = x_0;
f_new = f(:,end);
M = (L + R) / 2;
l = L;
r = R;
q = w(j);
s = a(j);
if (R(j) - L(j)) < 1.1*w(j)
  while 1
    M(j) = (l(j) + r(j))/2;
    
    [y_new, f_new] = getY(gp,xx,yy,z, f_new, M, eta, g);
    
    if s == 0 || y < y_new
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
  
  [y_new_ll] = getY(gp,xx,yy,z,f_new,tmp_ll, eta, g);  
  [y_new_rr] = getY(gp,xx,yy,z,f_new,tmp_rr, eta, g);
  
  if y >= y_new_ll
    ll(j) = ll(j) + q;
  end
  if y >= y_new_rr
    rr(j) = rr(j) - q;
  end
end % while
x_new(j) = ll(j) + rr(j) - x_0(j);
[y_new, f_new] = getY(gp,xx,yy,z,f_new, x_new, eta, g);
% y_new = -f(x_new,varargin{:});
if x_new(j) < l(j) || x_new(j) > r(j) || y >= y_new
  x_new(j) = x_0(j);
  rej = rej + 1;
  rej_old = 1;
else
  rej_old = 0;
  f(:,end+1) = f_new;
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function [x_new, f_new] = shrinkage(f,y,x_0,w,L,R,method,j,maxiter,um,xx,yy,gp,z,eta,g,varargin);
%function [x_new, y_new] = shrinkage(f,y,x_0,w,L,R,method,j,maxiter,um,varargin);
%    
% Shrinkage with acceptance-check for doubling scheme
% - acceptance-check is skipped if the distribution is defined
%   to be unimodal by the user

iter = 0;
x_new = x_0;
l = L(j);
r = R(j);

f_new = f(:,end);
% [y, tmp, eta, g] = getY(gp, xx, yy, z, f_new, x_0, []);

while 1
  x_new(j) = l + (r-l).*rand;
  
  [y_new, f_new] = getY(gp, xx, yy, z, [], x_new, eta, g);
  
  if strcmp(method,'doubling')
    if y < y_new && (um || accept(f,y,x_0,x_new,w,L,R,j,varargin{:}))
      break;
    end
  else
    if y < y_new
%       f(:,end+1) = f_new;
      break;
    end
  end % if strcmp
  if x_new(j) < x_0(j)
    l = x_new(j);
  else
    r = x_new(j);
  end % if
  if abs(l-r) < 1e-8
    error('bug')
  end
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

function [L,R] = stepping_out(f,y,x_0,L,R,w,m,j,mmin,mmax,di,um,xx,yy,gp,z,eta,g,varargin);
%function [L,R] = stepping_out(f,y,x_0,L,R,w,m,j,mmin,mmax,di,um,varargin);
%
% Stepping-out procedure

f_new = f(:,end);
x_new = x_0;
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
  y_new = getY(gp,xx,yy,z,f_new, L, eta, g);
  while y < y_new
%   while y < -f(L,varargin{:})
    L(j) = L(j) - w(j);
    if L(j) < mmin(j)
      L(j) = mmin(j);
      if di
        fprintf('Underflow! (L:%d)\n',j);
      end
      break;
    else
      y_new = getY(gp,xx,yy,z,f_new, L, eta, g);
    end
  end
  y_new = getY(gp,xx,yy,z,f_new, R, eta, g);
  while y < y_new
%   while y < -f(R,varargin{:})
    R(j) = R(j) + w(j);
    if R(j) > mmax(j)
      R(j) = mmax(j);
      if di
        fprintf('Overflow! (R:%d)\n',j);
      end
      break;
    else
        y_new = getY(gp,xx,yy,z,f_new, R, eta, g);
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
  y_new = getY(gp,xx,yy,z,f_new, L, eta, g);
  while J > 0 && y < y_new
%   while J > 0 && y < -f(L,varargin{:})
    L(j) = L(j) - w(j);
    if L(j) < mmin(j)
      L(j) = mmin(j);
      if di
        fprintf('Underflow! (L:%d)\n',j);
      end
      break;
    end
    y_new = getY(gp,xx,yy,z,f_new, L, eta, g);
    J = J - 1;
  end
  y_new = getY(gp,xx,yy,z,f_new, R, eta, g);
  while K > 0 && y < y_new
%   while K > 0 && y < -f(R,varargin{:})
    R(j) = R(j) + w(j);
    if R(j) > mmax(j)
      R(j) = mmax(j);
      if di
        fprintf('Overflow! (R:%d)\n',j);
      end
      break;
    end
    y_new = getY(gp,xx,yy,z,f_new, R, eta, g);
    K = K - 1;
  end
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

function [L,R] = doubling(f,y,x_0,L,R,w,m,j,mmin,mmax,di,um,xx,yy,gp,z,eta,g,varargin)
%function [L,R] = doubling(f,y,x_0,L,R,w,p,j,mmin,mmax,di,um,varargin);
%
% Doubling scheme for slice sampling
f_new = f(:,end);
x_new = x_0;
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
  AL = getY(gp,xx,yy,z,f_new, L, eta, g);
  AR = getY(gp,xx,yy,z,f_new, R, eta, g);
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
      AL = getY(gp,xx,yy,z,f_new, L, eta, g);
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
      AR = getY(gp,xx,yy,z,f_new, R, eta, g);
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
  AL = getY(gp,xx,yy,z,f_new, L, eta, g);
  AR = getY(gp,xx,yy,z,f_new, R, eta, g);
  while K > 0 && (y < AL || y < AR)
    if rand < 1/2
      L(j) = L(j) - (R(j)-L(j));
      if L(j) < mmin(j)
        L(j) = mmin(j);
        if di
          fprintf('Underflow! (L:%d)\n',j);
        end
      end
      AL = getY(gp,xx,yy,z,f_new, L, eta, g);
    else
      R(j) = R(j) + (R(j)-L(j));
      if R(j) > mmax(j)
        R(j) = mmax(j);
        if di
          fprintf('Overflow! (R:%d)\n',j);
        end
      end
      AR = getY(gp,xx,yy,z,f_new, R, eta, g);
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

function [y, f_new, eta, g] = getY(gp, xx, yy, z, f, w, eta, g)

if isempty(f) && (isempty(eta) || isempty(g))
  error('Must provide either current latent values f to get treshold or eta & g to get new latent values')
end

gp = gp_unpak(gp, w);
if ~isfield(gp.lik, 'nondiagW') || ismember(gp.lik.type, {'LGP' 'LGPC'}) 
  [K, C] = gp_trcov(gp, xx);
else
  if ~isfield(gp.lik,'xtime')
    nl=[0 repmat(size(yy,1), 1, length(gp.comp_cf))];
  else
    xtime=gp.lik.xtime;
    nl=[0 size(gp.lik.xtime,1) size(yy,1)];
  end
  nl=cumsum(nl);
  nlp=length(nl)-1;
  
  K = zeros(nl(end));
  for i1=1:nlp
    if i1==1 && isfield(gp.lik, 'xtime')
      K((1+nl(i1)):nl(i1+1),(1+nl(i1)):nl(i1+1)) = gp_trcov(gp, xtime, gp.comp_cf{i1});
    else
      K((1+nl(i1)):nl(i1+1),(1+nl(i1)):nl(i1+1)) = gp_trcov(gp, xx, gp.comp_cf{i1});
    end
  end
  C=K;
end
% for ii=1:size(yy,1)
%   [tmp,tmp, m2(ii,:)] = gp.lik.fh.tiltedMoments(gp.lik, yy, ii, C(ii,ii), 0, z);
% end
% S = diag(1./(1./m2 - 1./diag(C)));
S = 10*eye(size(K));
if isempty(eta) || isempty(g)
  g = mvnrnd(f,S)';
end
R = S-S*((S+K)\S);
R = (R+R')./2;
LR = chol(R,'lower');
m = R*(S\g);
if isempty(eta) || isempty(g)
  eta = LR\(f-m);
  f_new = [];
  tr = 1;   % return treshold
else
  f_new = LR*eta + m;
  tr = 0;   % return y for treshold comparison
end

% Log prior for proposed hyperparameters
lp = 0;
for i3=1:length(gp.cf)
  gpcf = gp.cf{i3};
  lp = lp + gpcf.fh.lp(gpcf);
end
if isfield(gp, 'lik') && isfield(gp.lik, 'p')
  likelih = gp.lik;
  lp = lp + likelih.fh.lp(likelih);
end

if tr
  % return treshold
  y = log(rand(1)) + gp.lik.fh.ll(gp.lik, yy, f, z) + mnorm_lpdf (g', 0, C + S) + lp;
else
  % return comparison value with proposed parameters
  y = gp.lik.fh.ll(gp.lik, yy, f_new, z) + mnorm_lpdf (g', 0, C + S) + lp;
end


function opt = ssls_opt(opt)
% Default opt for surrogate sls.
% fsamples - number of latent samples per hyperparameter sample
  
if ~isfield(opt, 'fsamples')
  opt.fsamples = 2;
end
if nargin < 1
  opt=[];
end
if nargin < 1
  opt=[];
end

if ~isfield(opt,'nsamples')
  opt.nsamples = 1;
end
if ~isfield(opt,'nomit')
  opt.nomit = 0;
end
if ~isfield(opt,'display')
  opt.display = 0;
end
if ~isfield(opt,'method')
  opt.method = 'multi';
end
if ~isfield(opt,'overrelaxation')
  opt.overrelaxation = 0;
elseif opt.overrelaxation == 1 && (strcmp(opt.method,'doubling') || strcmp(opt.method,'minmax'))
  opt.method = 'stepping';
end
if ~isfield(opt,'alimit')
  opt.alimit = 4;
end
if ~isfield(opt,'wsize')
  opt.wsize = 2;
end
if ~isfield(opt,'mlimit')
  opt.mlimit = 4;
end
if ~isfield(opt,'maxiter')
  opt.maxiter = 50;
end
if ~isfield(opt,'plimit')
  opt.plimit = 2;
end
if ~isfield(opt,'unimodal')
  opt.unimodal = 0;
end
if ~isfield(opt,'mmlimits')
  opt.mmlimits = [opt.wsize-(opt.wsize*opt.mlimit); opt.wsize+(opt.wsize*opt.mlimit)];
end

