function w = flattopwin (L, sym)
% Author: Paul Kienzle <pkienzle@users.sf.net> (2004)
% This program is granted to the public domain.
%
% flattopwin(L, [periodic|symmetric])
%
% Return the window f(w):
%
%   f(w) = 1 - 1.93 cos(2 pi w) + 1.29 cos(4 pi w)
%            - 0.388 cos(6 pi w) + 0.0322cos(8 pi w)
%
% where w = i/(L-1) for i=0:L-1 for a symmetric window, or
% w = i/L for i=0:L-1 for a periodic window.  The default
% is symmetric.  The returned window is normalized to a peak
% of 1 at w = 0.5.
%
% This window has low pass-band ripple, but high bandwidth.
%
% According to [1]:
%
%    The main use for the Flat Top window is for calibration, due
%    to its negligible amplitude errors.
%
% [1] Gade, S; Herlufsen, H; (1987) 'Use of weighting functions in DFT/FFT
% analysis (Part I)', Bruel & Kjaer Technical Review No.3.

  if nargin == 0 || nargin > 2
    help(mfilename);
  end % if

  divisor = L-1;
  if nargin > 1
    if strcmp(sym, 'periodic')
      divisor = L;
    elseif strcmp(sym, 'symmetric')
      divisor = L-1;
    else
      error('second argument must be ''periodic'' or ''symmetric''');
    end
  end % if
    
  x = 2*pi*(0:(L-1))'/divisor;
  if L==1
    w = 1;
  else
    % these coefficients come from wikipedia, and match better with what
    % matlab outputs (JM)
    a0 = 0.21557895;
    a1 = 0.41663158;
    a2 = 0.277263158;
    a3 = 0.083578947;
    a4 = 0.006947368;
    w = a0 - a1*cos(x) + a2*cos(2*x) - a3*cos(3*x) + ...
            a4*cos(4*x);
    %w = (1-1.93*cos(x)+1.29*cos(2*x)-0.388*cos(3*x)+0.0322*cos(4*x))/4.6402;
  end
end
