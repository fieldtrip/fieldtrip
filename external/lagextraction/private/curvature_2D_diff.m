function K = curvature_2D_diff(x,y)
%   CURVATURE_2D_DIFF   Compute curvature of a 2D curve
%       [K] = CURVATURE_2D_DIFF(X,Y)
%
%   Created by Alexandre Gramfort on 2008-12-16.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.

% $Id: curvature_2D_diff.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

% Smooth x and y with fix borders
niter = 5;
begx = x(1); endx = x(end);
begy = y(1); endy = y(end);
for ii=1:niter
    x = filter2([1/4;1/2;1/4],x(:),'valid')';
    y = filter2([1/4;1/2;1/4],y(:),'valid')';
    x = [begx,x,endx];
    y = [begy,y,endy];
end

% keyboard

% compute curvature
% dx = filter2([1/2;-1/2],x(:),'same')';
% dy = filter2([1/2;-1/2],y(:),'same')';
% 
% ddfilt = [1;-2;1];
% % ddfilt = [-1;16;-30;16;-1];
% ddfilt = ddfilt ./ sum(abs(ddfilt));
% 
% ddx = filter2(ddfilt,x(:),'same')';
% ddy = filter2(ddfilt,y(:),'same')';
% 
% s = sqrt( dx.^2 + dy.^2 );
% 
% K = abs(dx .* ddy - dy .* ddx) ./ s.^3;
% to_nan = 10;
% K(1:to_nan) = NaN;
% K(end-to_nan:end) = NaN
% return

% compute curvature
dx = compute_diff(x);
dy = compute_diff(y);

s = sqrt( dx.^2 + dy.^2 );
dx = dx./s; dy = dy./s; 
% compute dot product
p = length(dx);
K = zeros(p,1);
for k=1:p
    if k>1 && k<p
        K(k) = dx(k-1)*dx(k+1)+dy(k-1)*dy(k+1);
    elseif k==1
        K(k) = dx(k)*dx(k+1)+dy(k)*dy(k+1);
    elseif k==p
        K(k) = dx(k-1)*dx(k)+dy(k-1)*dy(k);
    end
end

function dx = compute_diff(x)

% compute_diff - compute central derivative of a vector.
%
% dx = compute_diff(x);

% central differences
h = 1; % step size
D1 = [x(2:end),x(end)];
D2 = [x(1),x(1:end-1)];
dx = (D1-D2)/(2*h); 
% y(1)=y(1)*2; y(end) = y(end)*2;
dx(1) = ( 4*x(2) - 3*x(1) - x(3) )/(2*h);
dx(end) = -( 4*x(end-1) - 3*x(end) - x(end-2) )/(2*h);

end %  function

end %  function