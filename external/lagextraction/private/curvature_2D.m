function K = curvature_2D(X,Y,degree)
% CURVATURE_2D estimate curvature of a curve defined by points
%
%   K = CURVATURE_2D(X, Y, DEGREE)
%
%   DEGREE is the degree of the polynom used for the interpolation
%

%
%   Created by Alexandre Gramfort on 2008-12-15.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.
%

% $Id: curvature_2D.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

if nargin<3
    % default values
    degree = 5; % degree of polynum used for interpolation
end

% Find a parametrization of the curve
t = cumsum(sqrt(diff(X).^2 + diff(Y).^2));
t = [0,t];
t = t / t(end); % normalize between 0 and 1

% compute coefficients of interpolation functions
x0 = polyfit(t, X, degree);
y0 = polyfit(t, Y, degree);

% compute coefficients of first and second derivatives. In the case of a
% polynom, it is possible to compute coefficient of derivative by
% multiplying with a matrix.
derive = diag(degree:-1:0);
xp = circshift(x0*derive, [0 1]);
yp = circshift(y0*derive, [0 1]);
xs = circshift(xp*derive, [0 1]);
ys = circshift(yp*derive, [0 1]);

% compute values of first and second derivatives for needed points 
xprime = polyval(xp, t);
yprime = polyval(yp, t);
xsec = polyval(xs, t);
ysec = polyval(ys, t);

% compute value of curvature
K = (xprime.*ysec - xsec.*yprime)./ ...
    power(xprime.*xprime + yprime.*yprime, 3/2);
