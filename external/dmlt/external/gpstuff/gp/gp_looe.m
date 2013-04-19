function eloo = gp_looe(w, gp, x, y)
%GP_LOOE  Evaluate the mean negative log leave-one-out predictive 
%         density, assuming Gaussian observation model
%
%  Description
%    LOOE = GP_LOOE(W, GP, X, Y) takes a parameter vector W,
%    Gaussian process structure GP, a matrix X of input vectors and
%    a matrix Y of targets, and evaluates the mean negative log
%    leave-one-out predictive density
%       LOOE  = - 1/n sum log p(Y_i | X, Y_{\i}, th)
%    where th represents the parameters (lengthScale,
%    magnSigma2...), X is inputs and Y is observations.
%
%  References:
%    S. Sundararajan and S. S. Keerthi (2001). Predictive
%    Approaches for Choosing Hyperparameters in Gaussian Processes. 
%    Neural Computation 13:1103-1118.
%
%  See also
%   GP_LOOG, GP_SET, GP_PAK, GP_UNPAK
%

% Copyright (c) 2008-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% Nothing to parse, but check the arguments anyway
ip=inputParser;
ip.FunctionName = 'GP_LOOE';
ip.addRequired('w', @(x) isvector(x) && isreal(x) && all(isfinite(x)));
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.parse(w, gp, x, y);

if isfield(gp,'mean') & ~isempty(gp.mean.meanFuncs)
  error('GP_LOOE: Mean functions not yet supported');
end

gp=gp_unpak(gp, w);
n = size(x,1);

% For single Gaussian process with Gaussian likelihood LOO
% log predictive density can be computed analytically.
% S. Sundararajan and S. S. Keerthi (2001). Predictive Approaches
% for Choosing Hyperparameters in Gaussian Processes. Neural
% Computation 13:1103-1118.

% gp_looprep returns b=C\y and iCv=diag(inv(C))
% using efficient computation for CS, FIC, PIC, and CS+FIC
[b,iCv]=gp_looprep(gp,x,y);
if isnan(b)
  eloo=NaN;
  return
end

% LOO-predictions
myy = y - b./iCv;
sigma2 = 1./iCv;
lpyt = (-0.5 * (log(2*pi) + log(sigma2) + (y-myy).^2./sigma2));
eloo = -sum(lpyt);

end
