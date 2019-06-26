function a=om_vangle(r,s)
% Returns an angle between two 3D vectors r,s in a range 0-pi

mag = sqrt(dot(r,r,2).*dot(s,s,2));
mag = max(mag,1e-30);
%if abs(mag)<1e-30, mag=1 ; end ; % avoid division by zero
a = acos(dot(r,s,2)./mag);
