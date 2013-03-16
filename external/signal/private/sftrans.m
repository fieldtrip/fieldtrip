% Copyright (C) 1999 Paul Kienzle
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.

% usage: [Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, W, stop)
%
% Transform band edges of a generic lowpass filter (cutoff at W=1)
% represented in splane zero-pole-gain form.  W is the edge of the
% target filter (or edges if band pass or band stop). Stop is true for
% high pass and band stop filters or false for low pass and band pass
% filters. Filter edges are specified in radians, from 0 to pi (the
% nyquist frequency).
%
% Theory: Given a low pass filter represented by poles and zeros in the
% splane, you can convert it to a low pass, high pass, band pass or
% band stop by transforming each of the poles and zeros individually.
% The following table summarizes the transformation:
%
% Transform         Zero at x                  Pole at x
% ----------------  -------------------------  ------------------------
% Low Pass          zero: Fc x/C               pole: Fc x/C
% S -> C S/Fc       gain: C/Fc                 gain: Fc/C
% ----------------  -------------------------  ------------------------
% High Pass         zero: Fc C/x               pole: Fc C/x
% S -> C Fc/S       pole: 0                    zero: 0
%                   gain: -x                   gain: -1/x
% ----------------  -------------------------  ------------------------
% Band Pass         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
%        S^2+FhFl   pole: 0                    zero: 0
% S -> C --------   gain: C/(Fh-Fl)            gain: (Fh-Fl)/C
%        S(Fh-Fl)   b=x/C (Fh-Fl)/2            b=x/C (Fh-Fl)/2
% ----------------  -------------------------  ------------------------
% Band Stop         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
%        S(Fh-Fl)   pole: ?sqrt(-FhFl)         zero: ?sqrt(-FhFl)
% S -> C --------   gain: -x                   gain: -1/x
%        S^2+FhFl   b=C/x (Fh-Fl)/2            b=C/x (Fh-Fl)/2
% ----------------  -------------------------  ------------------------
% Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
%      2 z-1        pole: -1                   zero: -1
% S -> - ---        gain: (2-xT)/T             gain: (2-xT)/T
%      T z+1
% ----------------  -------------------------  ------------------------
%
% where C is the cutoff frequency of the initial lowpass filter, Fc is
% the edge of the target low/high pass filter and [Fl,Fh] are the edges
% of the target band pass/stop filter.  With abundant tedious algebra,
% you can derive the above formulae yourself by substituting the
% transform for S into H(S)=S-x for a zero at x or H(S)=1/(S-x) for a
% pole at x, and converting the result into the form:
%
%    H(S)=g prod(S-Xi)/prod(S-Xj)
%
% The transforms are from the references.  The actual pole-zero-gain
% changes I derived myself.
%
% Please note that a pole and a zero at the same place exactly cancel.
% This is significant for High Pass, Band Pass and Band Stop filters
% which create numerous extra poles and zeros, most of which cancel.
% Those which do not cancel have a 'fill-in' effect, extending the
% shorter of the sets to have the same number of as the longer of the
% sets of poles and zeros (or at least split the difference in the case
% of the band pass filter).  There may be other opportunistic
% cancellations but I will not check for them.
%
% Also note that any pole on the unit circle or beyond will result in
% an unstable filter.  Because of cancellation, this will only happen
% if the number of poles is smaller than the number of zeros and the
% filter is high pass or band pass.  The analytic design methods all
% yield more poles than zeros, so this will not be a problem.
%
% References:
%
% Proakis & Manolakis (1992). Digital Signal Processing. New York:
% Macmillan Publishing Company.

% Author: Paul Kienzle <pkienzle@users.sf.net>

% 2000-03-01 pkienzle@kienzle.powernet.co.uk
%       leave transformed Sg as a complex value since cheby2 blows up
%       otherwise (but only for odd-order low-pass filters).  bilinear
%       will return Zg as real, so there is no visible change to the
%       user of the IIR filter design functions.
% 2001-03-09 pkienzle@kienzle.powernet.co.uk
%       return real Sg; don't know what to do for imaginary filters
function [Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, W, stop)

if (nargin ~= 5)
  usage('[Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, W, stop)');
end;

C = 1;
p = length(Sp);
z = length(Sz);
if z > p || p == 0
  error('sftrans: must have at least as many poles as zeros in s-plane');
end

if length(W)==2
  Fl = W(1);
  Fh = W(2);
  if stop
    % ----------------  -------------------------  ------------------------
    % Band Stop         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
    %        S(Fh-Fl)   pole: ?sqrt(-FhFl)         zero: ?sqrt(-FhFl)
    % S -> C --------   gain: -x                   gain: -1/x
    %        S^2+FhFl   b=C/x (Fh-Fl)/2            b=C/x (Fh-Fl)/2
    % ----------------  -------------------------  ------------------------
    if (isempty(Sz))
      Sg = Sg * real (1./ prod(-Sp));
    elseif (isempty(Sp))
      Sg = Sg * real(prod(-Sz));
    else
      Sg = Sg * real(prod(-Sz)/prod(-Sp));
    end
    b = (C*(Fh-Fl)/2)./Sp;
    Sp = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
    extend = [sqrt(-Fh*Fl), -sqrt(-Fh*Fl)];
    if isempty(Sz)
      Sz = [extend(1+rem([1:2*p],2))];
    else
      b = (C*(Fh-Fl)/2)./Sz;
      Sz = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
      if (p > z)
        Sz = [Sz, extend(1+rem([1:2*(p-z)],2))];
      end
    end
  else

    % ----------------  -------------------------  ------------------------
    % Band Pass         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
    %        S^2+FhFl   pole: 0                    zero: 0
    % S -> C --------   gain: C/(Fh-Fl)            gain: (Fh-Fl)/C
    %        S(Fh-Fl)   b=x/C (Fh-Fl)/2            b=x/C (Fh-Fl)/2
    % ----------------  -------------------------  ------------------------
    Sg = Sg * (C/(Fh-Fl))^(z-p);
    b = Sp*((Fh-Fl)/(2*C));
    Sp = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
    if isempty(Sz)
      Sz = zeros(1,p);
    else
      b = Sz*((Fh-Fl)/(2*C));
      Sz = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
      if (p>z)
        Sz = [Sz, zeros(1, (p-z))];
      end
    end
  end
else
  Fc = W;
  if stop
    % ----------------  -------------------------  ------------------------
    % High Pass         zero: Fc C/x               pole: Fc C/x
    % S -> C Fc/S       pole: 0                    zero: 0
    %                   gain: -x                   gain: -1/x
    % ----------------  -------------------------  ------------------------
    if (isempty(Sz))
      Sg = Sg * real (1./ prod(-Sp));
    elseif (isempty(Sp))
      Sg = Sg * real(prod(-Sz));
    else
      Sg = Sg * real(prod(-Sz)/prod(-Sp));
    end
    Sp = C * Fc ./ Sp;
    if isempty(Sz)
      Sz = zeros(1,p);
    else
      Sz = [C * Fc ./ Sz];
      if (p > z)
        Sz = [Sz, zeros(1,p-z)];
      end
    end
  else
    % ----------------  -------------------------  ------------------------
    % Low Pass          zero: Fc x/C               pole: Fc x/C
    % S -> C S/Fc       gain: C/Fc                 gain: Fc/C
    % ----------------  -------------------------  ------------------------
    Sg = Sg * (C/Fc)^(z-p);
    Sp = Fc * Sp / C;
    Sz = Fc * Sz / C;
  end
end
