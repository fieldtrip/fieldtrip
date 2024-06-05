function [A,D] = spm_matrix(P, order)
% Return an affine transformation matrix
% FORMAT [A] = spm_matrix(P [,order])
% P(1)  - x translation
% P(2)  - y translation
% P(3)  - z translation
% P(4)  - x rotation about - {pitch} (radians)
% P(5)  - y rotation about - {roll}  (radians)
% P(6)  - z rotation about - {yaw}   (radians)
% P(7)  - x scaling
% P(8)  - y scaling
% P(9)  - z scaling
% P(10) - x affine
% P(11) - y affine
% P(12) - z affine
%
% order - application order of transformations [Default: 'T*R*Z*S']
%
% A     - affine transformation matrix
%__________________________________________________________________________
%
% spm_matrix returns a matrix defining an orthogonal linear (translation,
% rotation, scaling or affine) transformation given a vector of
% parameters (P).  By default, the transformations are applied in the
% following order (i.e., the opposite to which they are specified):
%
% 1) shear
% 2) scale (zoom)
% 3) rotation - yaw, roll & pitch
% 4) translation
%
% This order can be changed by calling spm_matrix with a string as a
% second argument. This string may contain any valid MATLAB expression
% that returns a 4x4 matrix after evaluation. The special characters 'S',
% 'Z', 'R', 'T' can be used to reference the transformations 1)-4)
% above. The default order is 'T*R*Z*S', as described above.
%
% SPM uses a PRE-multiplication format i.e. Y = A*X where X and Y are 4 x n
% matrices of n coordinates.
%__________________________________________________________________________
%
% See also: spm_imatrix.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1994-2022 Wellcome Centre for Human Neuroimaging


%-Special case: translation only
%--------------------------------------------------------------------------
if numel(P) == 3
    A = eye(4);
    A(1:3,4) = P(:);
    return;
end
    
%-Pad P with 'null' parameters
%--------------------------------------------------------------------------
q  = [0 0 0 0 0 0 1 1 1 0 0 0];
P  = [P q((length(P) + 1):12)];

%-Translation / Rotation / Scale / Shear
%--------------------------------------------------------------------------
T  =   [1   0   0   P(1);
        0   1   0   P(2);
        0   0   1   P(3);
        0   0   0   1];

R1  =  [1   0           0           0;
        0   cos(P(4))   sin(P(4))   0;
        0  -sin(P(4))   cos(P(4))   0;
        0   0           0           1];

R2  =  [cos(P(5))   0   sin(P(5))   0;
        0           1   0           0;
       -sin(P(5))   0   cos(P(5))   0;
        0           0   0           1];

R3  =  [cos(P(6))   sin(P(6))   0   0;
       -sin(P(6))   cos(P(6))   0   0;
        0           0           1   0;
        0           0           0   1];

R   = R1*R2*R3;

Z   =  [P(7)   0       0       0;
        0      P(8)    0       0;
        0      0       P(9)    0;
        0      0       0       1];

S   =  [1      P(10)   P(11)   0;
        0      1       P(12)   0;
        0      0       1       0;
        0      0       0       1];

%-Affine transformation matrix
%--------------------------------------------------------------------------
if nargin < 2
    A = T*R*Z*S;
    if nargout>=2
        D = diff_matrix(P);
    end
else
    A = eval(sprintf('%s;', order));
    if ~isnumeric(A) || ~isequal(size(A),[4 4])
        error('Invalid order expression ''%s''.', order);
    end
end

function D = diff_matrix(p)
% Differentiate the matrix with respect to the parameters
n  = numel(p);
p0 = [0 0 0  0 0 0  1 1 1  0 0 0];
p0(1:n) = p;
p  = p0;

c4 = cos(p(4)); s4 = sin(p(4));
c5 = cos(p(5)); s5 = sin(p(5));
c6 = cos(p(6)); s6 = sin(p(6));
p7  = p(7);  p8  = p(8);  p9  = p(9);
p10 = p(10); p11 = p(11); p12 = p(12);

D = zeros(4,4,12);
D(1,4,1) = 1;
D(2,4,2) = 1;
D(3,4,3) = 1;

D(1:3,1:3,4) = [...
                    0,                                                 0,                                                                0
p7*(s4*s6 - c4*c6*s5), p7*p10*(s4*s6 - c4*c6*s5) - p8*(c6*s4 + c4*s5*s6), p9*c4*c5 + p7*p11*(s4*s6 - c4*c6*s5) - p8*p12*(c6*s4 + c4*s5*s6)
p7*(c4*s6 + c6*s4*s5), p7*p10*(c4*s6 + c6*s4*s5) - p8*(c4*c6 - s4*s5*s6), p7*p11*(c4*s6 + c6*s4*s5) - p9*c5*s4 - p8*p12*(c4*c6 - s4*s5*s6)];

D(1:3,1:3,5) = [...
   -p7*c6*s5,       - p8*s5*s6 - p7*p10*c6*s5,            p9*c5 - p8*p12*s5*s6 - p7*p11*c6*s5
-p7*c5*c6*s4, - p8*c5*s4*s6 - p7*p10*c5*c6*s4, - p9*s4*s5 - p7*p11*c5*c6*s4 - p8*p12*c5*s4*s6
-p7*c4*c5*c6, - p8*c4*c5*s6 - p7*p10*c4*c5*c6, - p9*c4*s5 - p7*p11*c4*c5*c6 - p8*p12*c4*c5*s6];

D(1:3,1:3,6) = [...
             -p7*c5*s6,                            p8*c5*c6 - p7*p10*c5*s6,                             p8*p12*c5*c6 - p7*p11*c5*s6
-p7*(c4*c6 - s4*s5*s6), -p8*(c4*s6 + c6*s4*s5) - p7*p10*(c4*c6 - s4*s5*s6), - p7*p11*(c4*c6 - s4*s5*s6) - p8*p12*(c4*s6 + c6*s4*s5)
 p7*(c6*s4 + c4*s5*s6),  p8*(s4*s6 - c4*c6*s5) + p7*p10*(c6*s4 + c4*s5*s6),   p7*p11*(c6*s4 + c4*s5*s6) + p8*p12*(s4*s6 - c4*c6*s5)];

D(1:3,1:3,7) = [...
            c5*c6,               p10*c5*c6,               p11*c5*c6
-c4*s6 - c6*s4*s5, -p10*(c4*s6 + c6*s4*s5), -p11*(c4*s6 + c6*s4*s5)
 s4*s6 - c4*c6*s5,  p10*(s4*s6 - c4*c6*s5),  p11*(s4*s6 - c4*c6*s5)];

D(1:3,1:3,8) = [...
0,              c5*s6,               p12*c5*s6
0,   c4*c6 - s4*s5*s6,  p12*(c4*c6 - s4*s5*s6)
0, - c6*s4 - c4*s5*s6, -p12*(c6*s4 + c4*s5*s6)];

D(1:3,1:3,9) = [...
0, 0,    s5
0, 0, c5*s4
0, 0, c4*c5];

D(1:3,1:3,10) = [...
0,               p7*c5*c6, 0
0, -p7*(c4*s6 + c6*s4*s5), 0
0,  p7*(s4*s6 - c4*c6*s5), 0];

D(1:3,1:3,11) = [...
0, 0,               p7*c5*c6
0, 0, -p7*(c4*s6 + c6*s4*s5)
0, 0,  p7*(s4*s6 - c4*c6*s5)];

D(1:3,1:3,12) = [...
0, 0,               p8*c5*s6
0, 0,  p8*(c4*c6 - s4*s5*s6)
0, 0, -p8*(c6*s4 + c4*s5*s6)];

D = D(:,:,1:n);



