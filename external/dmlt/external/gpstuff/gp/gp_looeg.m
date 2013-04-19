function [e, g] = gp_looeg(w, gp, x, y, varargin)
%GP_LOOEG  Evaluate the mean negative log leave-one-out predictive density
%
%  Description
%    [E, G] = GP_LOOEG(W, GP, X, Y, OPTIONS) takes a Gaussian
%    process structure GP together with a matrix X of input vectors
%    and a matrix Y of targets, and evaluates the energy function E
%    and its gradient G (only with Gaussian likelihood). Each row
%    of X corresponds to one input vector and each row of Y
%    corresponds to one target vector.
%
%    The energy is the mean negative log
%    leave-one-out predictive density
%       LOOE  = - 1/n sum log p(Y_i | X, Y_{\i}, th)
%    where th represents the parameters (lengthScale,
%    magnSigma2...), X is inputs and Y is observations.
%
%    For non-Gaussian models EP leave-one-out is used for energy
%    (see GPEP_LOOE), but no gradients are yet implemented, and
%    thus gradient-free optimization function has to be used.
%
%    OPTIONS is optional parameter-value pair
%      z - optional observed quantity in triplet (x_i,y_i,z_i)
%          Some likelihoods may use this. For example, in case of
%          Poisson likelihood we have z_i=E_i, that is, expected
%          value for ith case.
%
%  See also
%    GP_LOOE, GP_LOOG, GPEP_LOOE
%

% Copyright (c) 2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% Nothing to parse, but check the arguments anyway
ip=inputParser;
ip.FunctionName = 'GP_LOOEG';
ip.addRequired('w', @(x) isvector(x) && isreal(x) && all(isfinite(x)));
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.parse(w, gp, x, y);

if isfield(gp,'mean') & ~isempty(gp.mean.meanFuncs)
  error('GP_LOOEG: Mean functions not yet supported');
end

gp=gp_unpak(gp, w);
n = size(x,1);

% Single function for some optimization routines, no need for mydeal...
if isfield(gp.lik.fh,'trcov')
  % Gaussian likelihood
  if nargout<2
    % gp_looprep returns b=C\y and iCv=diag(inv(C))
    % using efficient computation for CS, FIC, PIC, and CS+FIC
    [b,iCv]=gp_looprep(gp,x,y);
    if isnan(b)
      e=NaN;
      return
    end
    myy = y - b./iCv;
    sigma2 = 1./iCv;
    lpyt = (-0.5 * (log(2*pi) + log(sigma2) + (y-myy).^2./sigma2));
    e = -sum(lpyt);
  else
    % gp_looprep returns b=C\y and iCv=diag(inv(C))
    % using efficient computation for CS, FIC, PIC, and CS+FIC
    % to avoid recomputation of these terms, both error and gradients are
    % computed below, instead of calling gp_looe and gp_loog
    [b,iCv,iC]=gp_looprep(gp,x,y);
    if isnan(b)
      e=NaN;
      g=NaN;
      return
    end
    myy = y - b./iCv;
    sigma2 = 1./iCv;
    lpyt = (-0.5 * (log(2*pi) + log(sigma2) + (y-myy).^2./sigma2));
    e = -sum(lpyt);
    
    % Get the gradients of the covariance matrices and gprior
    % from gpcf_* structures and evaluate the gradients
    ncf = length(gp.cf);
    g = [];
    i1=0;
    if ~isempty(strfind(gp.infer_params, 'covariance'))
      for i=1:ncf
        gpcf = gp.cf{i};
        gpcf.GPtype = gp.type;
        DKff = gpcf.fh.cfg(gpcf, x);
        
        % Evaluate the gradient with respect to covariance function parameters
        for i2 = 1:length(DKff)
          i1 = i1+1;  
          Z = iC*DKff{i2};
          Zb = Z*b;            
          g(i1) = - sum( (b.*Zb - 0.5*(1 + b.^2./iCv).*diag(Z*iC))./iCv );
        end
      end
    end

    % Evaluate the gradient from Gaussian likelihood function
    if ~isempty(strfind(gp.infer_params, 'likelihood'))
      if isfield(gp.lik.fh,'trcov')
        DCff = gp.lik.fh.cfg(gp.lik, x);
        for i2 = 1:length(DCff)
          i1 = i1+1;
          Z = iC*eye(n,n).*DCff{i2};
          Zb = Z*b;            
          g(i1) = - sum( (b.*Zb - 0.5*(1 + b.^2./iCv).*diag(Z*iC))./iCv );
        end
      end
    end
    
  end
else
  % non-Gaussian likelihood
  e=gp.fh.looe(w, gp, x, y, varargin{:});
  if nargout>1
    g=gp.fh.loog(w,gp,x,y);
  end
  if isnan(e)
    e=realmax;
  end
end
