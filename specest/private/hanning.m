function [tap] = hanning(n, str)

%HANNING   Hanning window.
%   HANNING(N) returns the N-point symmetric Hanning window in a column
%   vector.  Note that the first and last zero-weighted window samples
%   are not included.
%
%   HANNING(N,'symmetric') returns the same result as HANNING(N).
%
%   HANNING(N,'periodic') returns the N-point periodic Hanning window,
%   and includes the first zero-weighted window sample.
%
%   NOTE: Use the HANN function to get a Hanning window which has the 
%          first and last zero-weighted samples. 
%
%   See also BARTLETT, BLACKMAN, BOXCAR, CHEBWIN, HAMMING, HANN, KAISER
%   and TRIANG.
%   
%   This is a drop-in replacement to bypass the signal processing toolbox

% Copyright (c) 2010, Jan-Mathijs Schoffelen, DCCN Nijmegen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin==1,
  str = 'symmetric';
end
  
switch str,
case 'periodic'
   % Includes the first zero sample
   tap = [0; hanningX(n-1)];
case 'symmetric'
   % Does not include the first and last zero sample
   tap = hanningX(n);
end

function tap = hanningX(n)

% compute taper
N   = n+1;
tap = 0.5*(1-cos((2*pi*(1:n))./N))';

% make symmetric
halfn = floor(n/2);
tap( (n+1-halfn):n ) = flipud(tap(1:halfn));




