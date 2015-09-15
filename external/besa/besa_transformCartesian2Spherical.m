function [theta, phi, r] = besa_transformCartesian2Spherical(X, Y, Z)
% BESA_TRANSFORMCARTESIAN2SPHERICAL takes the cartesian coordinates X, Y 
% and Z of a 3D point and transforms them to the spherical coordinates 
% theta, phi and r. 
%
% Parameters:
%     [X]
%         The X-coordinate of the current 3D point. It should point to the
%         right.
% 
%     [Y]
%         The Y-coordinate of the current 3D point. It should point
%         foreward.
% 
%     [Z]
%         The Z-coordinate of the current 3D point. It should point up.
% 
% 
% Return:
%     [theta] 
%         The azimuth angle with the vertical z-axis (0 degree) in the 
%         x-z-plane, +90 degree lie in the positive x-axis.
% 
%     [phi] 
%         The latitude angle in the horizontal x-y-plane 
%         (counter-clockwise).
%
%     [r] 
%         The radius.

% Copyright (C) 2015, BESA GmbH
%
% File name: besa_transformCartesian2Spherical.m
%
% Author: Todor Jordanov
% Created: 2015-07-29

SquareOfRadiusXY = X*X + Y*Y;
RadiusXY = sqrt(SquareOfRadiusXY);

if(RadiusXY == 0.0 && Z == 0.0)

    theta = 0.0;

else

    theta = atan2d(RadiusXY, Z);

end

if(X==0.0 && Y==0.0)

    phi = 0.0;

else

    phi = atan2d(Y, X);

end

% square of total radius
SquareOfRadiusXY = SquareOfRadiusXY + Z*Z;

r = sqrt(SquareOfRadiusXY);

% set phi & theta to BESA ranges
if(phi > 90.)

    phi = phi - 180.0;
    theta = -theta;

elseif(phi < -90.) 

    phi = phi+180.0;
    theta = -theta;

end
