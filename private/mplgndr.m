function plgndr = mplgndr ( degree, order, x )

% MPLGNDR associated Legendre functions
%
% y = mplgndr(n,k,x) computes the values of the associated Legendre
% functions of order K up to degree N.
% 
% The input x can be a matrix, and the result is of size numel(x) by N+1.
% The i-th column is the associated Legendre function of order K and
% degree i-1.

% Copyright (C) 2002, Robert Oostenveld
% Copyright (C) 2021, Ricardo Bruña
% 
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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


% Checks the input.
if order < 0 || round ( order ) ~= order
  error ( 'The order must be a positive integer or zero.' )
end
if degree < order || round ( degree ) ~= degree
  error ( 'The degree must be a positive integer equal or greated than m.' )
end
if any ( abs ( x ) > 1 )
  error ( 'abs (x) vannot be greater than 1.' )
end


% Initializes the output.
plgndr = zeros ( numel ( x ), degree + 1 );


% Gets the solution at degree = order analitically.
root   = sqrt ( 1 - x (:) .* x (:) );
plgndr ( :, order + 1 ) = prod ( 1: 2: 2 * order - 1 ) * root .^ order;

% Corrects the sign.
if rem ( abs ( order ), 2 )
  plgndr ( :, order + 1 ) = -plgndr ( :, order + 1 );
end

% If the degree equals the order, returns.
if degree == order, return, end


% Uses a simplified version of the formula for the second iteration.
plgndr ( :, order + 2 ) = x (:) .* ( 2 * ( order + 1 ) - 1 ) .* plgndr ( :, order + 1 );


% Uses the iterative formula for the rest of iterations.
for dindex = order + 2: degree
  
  % Calculates the value of the polynomial in the current iteration.
  plgndr ( :, dindex + 1 ) = ( x (:) .* ( 2 * dindex - 1 ) .* plgndr ( :, dindex ) - ( dindex + order - 1 ) .* plgndr ( :, dindex - 1 ) ) ./ ( dindex - order );
end
