function P = spm_imatrix(M)
% returns the parameters for creating an affine transformation
% FORMAT P = spm_imatrix(M)
% M      - Affine transformation matrix
% P      - Parameters (see spm_matrix for definitions)
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Stefan Kiebel
% $Id: spm_imatrix.m 1143 2008-02-07 19:33:33Z spm $


% Translations and zooms
%-----------------------------------------------------------------------
R         = M(1:3,1:3);
C         = chol(R'*R);
P         = [M(1:3,4)' 0 0 0  diag(C)'  0 0 0];
if det(R)<0, P(7)=-P(7);end % Fix for -ve determinants

% Shears
%-----------------------------------------------------------------------
C         = diag(diag(C))\C;
P(10:12)  = C([4 7 8]);
R0        = spm_matrix([0 0 0  0 0 0 P(7:12)]);
R0        = R0(1:3,1:3);
R1        = R/R0;

% This just leaves rotations in matrix R1
%-----------------------------------------------------------------------
%[          c5*c6,           c5*s6, s5]
%[-s4*s5*c6-c4*s6, -s4*s5*s6+c4*c6, s4*c5]
%[-c4*s5*c6+s4*s6, -c4*s5*s6-s4*c6, c4*c5]

P(5) = asin(rang(R1(1,3)));
if (abs(P(5))-pi/2)^2 < 1e-9,
    P(4) = 0;
    P(6) = atan2(-rang(R1(2,1)), rang(-R1(3,1)/R1(1,3)));
else
    c    = cos(P(5));
    P(4) = atan2(rang(R1(2,3)/c), rang(R1(3,3)/c));
    P(6) = atan2(rang(R1(1,2)/c), rang(R1(1,1)/c));
end;
return;

% There may be slight rounding errors making b>1 or b<-1.
function a = rang(b)
a = min(max(b, -1), 1);
return;
