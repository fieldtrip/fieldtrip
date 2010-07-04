function a = angle_dist_nc(p2, p1)
%function to calculate the angle on the sphere going from from OP1 to OP2
%for a companion function, see RotateAboutAxis
%equivalent in SUMA code is macro : SUMA_ANGLE_DIST_NC
m_cr = cross(p1, p2);
a = atan2(sqrt(sum(m_cr.^2)),dot(p2, p1)); 

return;
