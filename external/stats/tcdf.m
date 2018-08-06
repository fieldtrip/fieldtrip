function p = tcdf(t,v,uflag)

% TCDF Student's T cumulative distribution function (cdf).
%
% P = TCDF(X,V) computes the cdf for Student's T distribution
% with V degrees of freedom, at the values in X. V must be a
% scalar or have the same size as T.
%
% P = TCDF(X,V,'upper') computes it for the upper tail instead
% of the lower tail.
%
% This is an alternative to the TCDF function that is implemented
% in the Matlab statistics toolbox. This version originates from
% http://www.statsci.org/matlab/statbox.html and originally was called TP.
% It has been renamed to TCDF for drop-in compatibility with the Matlab
% version.
%
% Gordon Smyth, University of Queensland, gks@maths.uq.edu.au
% 3 Apr 97
%
% NaN compatible - Markus Bauer and Eric Maris, FCDC
% 27 Jan 2005
%
% fixed bug concerning NaN compatibility
% 21 Aug 2006, Markus Siegel
%
% added support for upper tail, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3045
% 13 Jan 2016, Robert Oostenveld

if v <= 0, error('Degrees of freedom must be positive.'); end;

% resize v if necessary
if all(size(v)==1)
  v = ones(size(t))*v;
end;

%check for NaN's - don't do calculations on them, give those out as NaNs
if any( not(isfinite(t(:))) | not(isfinite(v(:))) )
  sel = find(isfinite(t) & isfinite(v));
  x=nan(size(t));
  p=nan(size(t));
  x(sel) = t(sel).^2 ./ (v(sel) + t(sel).^2) ;
  p(sel) = 0.5 .* ( 1 + sign(t(sel)) .* betainc( x(sel), 0.5, 0.5*v(sel) ) );
else
  x = t.^2 ./ (v + t.^2) ;
  p = 0.5 .* ( 1 + sign(t) .* betainc( x, 0.5, 0.5*v ) );
end;

if nargin==3
  if strcmp(uflag, 'upper')
    % compute it for the upper tail instead of the lower tail
    p = 1-p;
  else
    error('incorrect specification of uflag');
  end
end


