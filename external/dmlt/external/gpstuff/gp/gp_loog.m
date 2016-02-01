function gloo = gp_loog(w, gp, x, y, varargin)
%GP_LOOG  Evaluate the gradient of the mean negative log 
%         leave-one-out predictive density, assuming Gaussian 
%         observation model
%
%   Description
%     LOOG = GP_LOOG(W, GP, X, Y) takes a parameter vector W,
%     Gaussian process structure GP, a matrix X of input vectors
%     and a matrix Y of targets, and evaluates the gradient of the
%     mean negative log leave-one-out predictive density (see
%     GP_LOOE).
%
%   References:
%     S. Sundararajan and S. S. Keerthi (2001). Predictive
%     Approaches for Choosing Hyperparameters in Gaussian Processes. 
%     Neural Computation 13:1103-1118.
%
%  See also
%    GP_LOOE, GP_SET, GP_PAK, GP_UNPAK
%

% Copyright (c) 2008 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% Nothing to parse, but check the arguments anyway
ip=inputParser;
ip.FunctionName = 'GP_LOOG';
ip.addRequired('w', @(x) isvector(x) && isreal(x) && all(isfinite(x)));
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.parse(w, gp, x, y);

if isfield(gp,'mean') & ~isempty(gp.mean.meanFuncs)
  error('GP_LOOE: Mean functions not yet supported');
end

gp=gp_unpak(gp, w);
ncf = length(gp.cf);
n = size(x,1);

g = [];
gloo = [];

% For single Gaussian process with Gaussian likelihood LOO
% gradient of log predictive density can be computed analytically.
% S. Sundararajan and S. S. Keerthi (2001). Predictive Approaches
% for Choosing Hyperparameters in Gaussian Processes. Neural
% Computation 13:1103-1118.

% gp_looprep returns b=C\y and iCv=diag(inv(C))
% using efficient computation for CS, FIC, PIC, and CS+FIC
[b,iCv,iC]=gp_looprep(gp,x,y);
if isnan(b)
  gloo=NaN;
  return
end

% Get the gradients of the covariance matrices and gprior
% from gpcf_* structures and evaluate the gradients
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
      gloo(i1) = - sum( (b.*Zb - 0.5*(1 + b.^2./iCv).*diag(Z*iC))./iCv );
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
      gloo(i1) = - sum( (b.*Zb - 0.5*(1 + b.^2./iCv).*diag(Z*iC))./iCv );
    end
  end
end

end
