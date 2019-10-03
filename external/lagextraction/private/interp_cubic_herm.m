function y = interp_cubic_herm ( xi, fi, fpi, x )

% INTERP_CUBIC_HERM   evaluate the Hermite cubic interpolant associated with a 
%             given set of interpolating points, function values and 
%             derivative values at a specified set of values for the
%             independent variable
%
%     calling sequences:
%             y = interp_cubic_herm ( xi, fi, fpi, x )
%
%     inputs:
%             xi      vector containing the interpolating points
%             fi      vector containing function values
%                     the i-th entry in this vector is the function
%                     value associated with the i-th entry in the 'xi'
%                     vector
%             fpi     vector containing derivative values
%                     the i-th entry in this vector is the derivative
%                     value associated with the i-th entry in the 'xi'
%                     vector
%             x       value(s) of independent variable at which to
%                     evaluate the Hermite cubic interpolant
%                     - may be a scalar or a vector
%
%     output:
%             y       value(s) of the Hermite cubic interpolant
%

% $Id: interp_cubic_herm.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

if ~issorted(xi)
    [xi,I] = sort(xi);
    fi = fi(I);
    fpi = fpi(I);
end

if ~issorted(x)
    x = sort(x);
end

p = length ( xi );
npts = length ( x );

for i = 1 : npts
    j = max ( find ( xi(1:p-1) < x(i) ) );
    if ( length ( j ) == 0 )
       j = 1;
    end

    hj = ( xi(j+1) - xi(j) );
    temp = ( x(i) - xi(j) ) / hj;
    phi = ( 1 + 2*temp ) * ( 1 - temp )^2;
    psi1 = temp * ( 1 - temp )^2;
    psi2 = temp^2 * ( 1 - temp );
    out(i) = fi(j+1) + phi * ( fi(j) - fi(j+1) ) + hj * ( psi1 * fpi(j) - psi2 * fpi(j+1) );
end

if ( nargout == 0 )
   disp ( out )
else
   y = out;
end