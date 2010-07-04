function [XYZr] = RotateAboutAxis(XYZ, nrm, phi)
%[XYZr] = RotateAboutAxis(XYZ, nrm, phi)
%yet another function to perform rotations about an axis
%XYZ: Matrix of XYZ coordinates, each row is an XYZ triplet
%nrm: axis of rotation
%phi: Clockwise rotation angle in radians
%XYZr: Rotated version of XYZ
%Note:
%Origin of coordinate system is assumed to be 0
%Function is equivalent to AxisRotate3D
%equation for rotation from http://mathworld.wolfram.com/RotationFormula.html
%see C macro in SUMA: SUMA_ROTATE_ABOUT_AXIS
% ZSS, Nov 05

% angle already in radians phi = phi./180.*pi;

nrm = nrm./sqrt(sum(nrm.^2));

XYZr = zeros(size(XYZ));
for (i=1:1:size(XYZr,1)),
   XYZr(i,:) = XYZ(i,:).*cos(phi) + nrm.*(dot(nrm,XYZ(i,:))).*(1-cos(phi)) - (cross(XYZ(i,:),nrm)).*sin(phi);
end

if (0),
   %one at a time
   for (i=1:1:size(XYZ,1)),
      dop = dot(nrm,XYZ(i,:));
      cro = cross(XYZ(i,:),nrm);
      cop = cos(phi);
      sip = sin(phi);
      XYZr2(i,1) = XYZ(i,1).*cop + nrm(1).*(dop).*(1-cop) - (cro(1)).*sip;
      XYZr2(i,2) = XYZ(i,2).*cop + nrm(2).*(dop).*(1-cop) - (cro(2)).*sip;
      XYZr2(i,3) = XYZ(i,3).*cop + nrm(3).*(dop).*(1-cop) - (cro(3)).*sip;
   end
   XYZr-XYZr2
end

return;

