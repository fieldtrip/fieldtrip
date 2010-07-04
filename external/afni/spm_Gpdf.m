function f = spm_Gpdf(x,h,l)
% Probability Density Function (PDF) of Gamma distribution
% FORMAT f = spm_Gpdf(g,h,l)
%
% x - Gamma-variate   (Gamma has range [0,Inf) )
% h - Shape parameter (h>0)
% l - Scale parameter (l>0)
% f - PDF of Gamma-distribution with shape & scale parameters h & l
%_______________________________________________________________________
%
% spm_Gpdf implements the Probability Density Function of the Gamma
% distribution.
%
% Definition:
%-----------------------------------------------------------------------
% The PDF of the Gamma distribution with shape parameter h and scale l
% is defined for h>0 & l>0 and for x in [0,Inf) by: (See Evans et al.,
% Ch18, but note that this reference uses the alternative
% parameterisation of the Gamma with scale parameter c=1/l)
%
%           l^h * x^(h-1) exp(-lx)
%    f(x) = ---------------------
%                gamma(h)
%
% Variate relationships: (Evans et al., Ch18 & Ch8)
%-----------------------------------------------------------------------
% For natural (strictly +ve integer) shape h this is an Erlang distribution.
%
% The Standard Gamma distribution has a single parameter, the shape h.
% The scale taken as l=1.
%
% The Chi-squared distribution with v degrees of freedom is equivalent
% to the Gamma distribution with scale parameter 1/2 and shape parameter v/2.
%
% Algorithm:
%-----------------------------------------------------------------------
% Direct computation using logs to avoid roundoff errors.
%
% References:
%-----------------------------------------------------------------------
% Evans M, Hastings N, Peacock B (1993)
%       "Statistical Distributions"
%        2nd Ed. Wiley, New York
%
% Abramowitz M, Stegun IA, (1964)
%       "Handbook of Mathematical Functions"
%        US Government Printing Office
%
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
%       "Numerical Recipes in C"
%        Cambridge
%_______________________________________________________________________
% @(#)spm_Gpdf.m	2.2 Andrew Holmes 99/04/26

%-Format arguments, note & check sizes
%-----------------------------------------------------------------------
if nargin<3, error('Insufficient arguments'), end

ad = [ndims(x);ndims(h);ndims(l)];
rd = max(ad);
as = [	[size(x),ones(1,rd-ad(1))];...
	[size(h),ones(1,rd-ad(2))];...
	[size(l),ones(1,rd-ad(3))]     ];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 & any(any(diff(as(xa,:)),1))
	error('non-scalar args must match in size'), end

%-Computation
%-----------------------------------------------------------------------
%-Initialise result to zeros
f = zeros(rs);

%-Only defined for strictly positive h & l. Return NaN if undefined.
md = ( ones(size(x))  &  h>0  &  l>0 );
if any(~md(:)), f(~md) = NaN;
	warning('Returning NaN for out of range arguments'), end

%-Degenerate cases at x==0: h<1 => f=Inf; h==1 => f=l; h>1 => f=0
ml = ( md  &  x==0  &  h<1 );
f(ml) = Inf;
ml = ( md  &  x==0  &  h==1 ); if xa(3), mll=ml; else mll=1; end
f(ml) = l(mll);

%-Compute where defined and x>0
Q  = find( md  &  x>0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qh=Q; else Qh=1; end
if xa(3), Ql=Q; else Ql=1; end

%-Compute
f(Q) = exp( (h(Qh)-1).*log(x(Qx)) +h(Qh).*log(l(Ql)) - l(Ql).*x(Qx)...
		-gammaln(h(Qh)) );
