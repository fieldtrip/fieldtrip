function [cc] = circumcenter(coil1,coil2,coil3)
 
% CIRCUMCENTER determines the position and orientation of the circumcenter
% of the three fiducial markers (MEG headposition coils). 
%
% Input: X,y,z-coordinates of the 3 coils [3 X N],[3 X N],[3 X N] where N
% is timesamples/trials.
%
% Output: X,y,z-coordinates of the circumcenter [1-3 X N], and the 
% orientations to the x,y,z-axes [4-6 X N].
%
% A. Stolk, 2012
 
% number of timesamples/trials
N = size(coil1,2);
 
%% x-, y-, and z-coordinates of the circumcenter
% use coordinates relative to point `a' of the triangle
xba = coil2(1,:) - coil1(1,:);
yba = coil2(2,:) - coil1(2,:);
zba = coil2(3,:) - coil1(3,:);
xca = coil3(1,:) - coil1(1,:);
yca = coil3(2,:) - coil1(2,:);
zca = coil3(3,:) - coil1(3,:);
 
% squares of lengths of the edges incident to `a'
balength = xba .* xba + yba .* yba + zba .* zba;
calength = xca .* xca + yca .* yca + zca .* zca;
 
% cross product of these edges
xcrossbc = yba .* zca - yca .* zba;
ycrossbc = zba .* xca - zca .* xba;
zcrossbc = xba .* yca - xca .* yba;
 
% calculate the denominator of the formulae
denominator = 0.5 ./ (xcrossbc .* xcrossbc + ycrossbc .* ycrossbc + zcrossbc .* zcrossbc);
 
% calculate offset (from `a') of circumcenter
xcirca = ((balength .* yca - calength .* yba) .* zcrossbc - (balength .* zca - calength .* zba) .* ycrossbc) .* denominator;
ycirca = ((balength .* zca - calength .* zba) .* xcrossbc - (balength .* xca - calength .* xba) .* zcrossbc) .* denominator;
zcirca = ((balength .* xca - calength .* xba) .* ycrossbc - (balength .* yca - calength .* yba) .* xcrossbc) .* denominator;
 
cc(1,:) = xcirca + coil1(1,:);
cc(2,:) = ycirca + coil1(2,:);
cc(3,:) = zcirca + coil1(3,:);
 
%% orientation of the circumcenter with respect to the x-, y-, and z-axis
% coordinates
v = [cc(1,:)', cc(2,:)', cc(3,:)'];
vx = [zeros(1,N)', cc(2,:)', cc(3,:)']; % on the x-axis
vy = [cc(1,:)', zeros(1,N)', cc(3,:)']; % on the y-axis
vz = [cc(1,:)', cc(2,:)', zeros(1,N)']; % on the z-axis
 
for j = 1:N
  % find the angles of two vectors opposing the axes
  thetax(j) = acos(dot(v(j,:),vx(j,:))/(norm(v(j,:))*norm(vx(j,:))));
  thetay(j) = acos(dot(v(j,:),vy(j,:))/(norm(v(j,:))*norm(vy(j,:))));
  thetaz(j) = acos(dot(v(j,:),vz(j,:))/(norm(v(j,:))*norm(vz(j,:))));
 
  % convert to degrees
  cc(4,j) = (thetax(j) * (180/pi));
  cc(5,j) = (thetay(j) * (180/pi));
  cc(6,j) = (thetaz(j) * (180/pi));
end

