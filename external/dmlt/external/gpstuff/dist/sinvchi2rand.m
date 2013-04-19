function r = sinvchi2rand(nu, s2, M, N)
% SINVCHI2RAND  Random matrices from scaled inverse-chi distribution
%
%  R = SINVCHI2RAND(NU, S2)
%  R = SINVCHI2RAND(NU, S2, M, N)
%
%  Returns a randon number/matrix R from scaled inverse-chi square 
%  distribution. Nu is the degrees of freedom and S2 is the scale 
%  squared. Parametrisation is according to Gelman et. al. (2004).

% Copyright (c) 1998-2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 2
  error('Too few arguments');
end
if nargin==2
    [M,N]=size(s2);
else
    if numel(s2)>1 || numel(nu)>1
        error('Arguments M and N can only be used if nu and s2 are scalars');
    end
    if nargin < 4
        N=1;
    end
end
r=nu.*s2./chi2rnd(nu,M,N);
