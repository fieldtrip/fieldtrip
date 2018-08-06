function sixdof = hom2six(M)
% function sixdof = hom2six(M)
% 
% Convert homogenous transformation (4x4 matrix) to 6-dof representation
% as used in SPM8. That is, this function assumes that
% M = Trans([x;y;z]) * RotX(a) * RotY(b) * RotZ(c)
% and returns the 6 parameters [x,y,z,a,b,c] that generate M.
%
% Note that this representation is not unique in case b = +/- pi/2.

% (C) 2010 S. Klanke

%Rotation part is Rx(a)*Ry(b)*Rz(c)
%[           cb*cc,           cb*sc,              sb]
%[ -sa*sb*cc-ca*sc, -sa*sb*sc+ca*cc,           sa*cb]
%[ -ca*sb*cc+sa*sc, -ca*sb*sc-sa*cc,           ca*cb]

if M(1,3) == 1   % sb = 1, cb = 0
%[  0, 0, 1]
%[ -sa*cc-ca*sc, -sa*sc+ca*cc, 0]
%[ -ca*cc+sa*sc, -ca*sc-sa*cc, 0]
% not unique -> set a = 0 - > ca=1, sa=0
%[  0   0  1]
%[-sc  cc  0]
%[-cc -sc  0]
	a = 0;
	b = pi/2;
	c = atan2(-M(3,2),M(2,2));
elseif M(1,3) == -1  % sb = -1, cb = 0
%[  0, 0, -1]
%[ +sa*cc-ca*sc, +sa*sc+ca*cc, 0]
%[ +ca*cc+sa*sc, +ca*sc-sa*cc, 0]
% not unique -> set a = 0 - > ca=1, sa=0
%[  0   0  -1]
%[-sc  cc   0]
%[ cc  sc   0]
    a = 0;
	b = -pi/2;
	c = atan2(M(3,2),M(2,2));
else
   % cos(b)~=0
   c = atan2(M(1,2),M(1,1));
   b = atan2(M(1,3),norm(M(1,1:2)));
   a = atan2(M(2,3),M(3,3));
end

sixdof = [M(1:3,4)' a b c];
