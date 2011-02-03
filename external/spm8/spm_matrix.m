function [A] = spm_matrix(P, order)
% returns an affine transformation matrix
% FORMAT [A] = spm_matrix(P, order)
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
% order (optional) application order of transformations.
%
% A     - affine transformation matrix
%___________________________________________________________________________
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
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id$


% pad P with 'null' parameters
%---------------------------------------------------------------------------
q  = [0 0 0 0 0 0 1 1 1 0 0 0];
P  = [P q((length(P) + 1):12)];

% default multiplication order if not specified
%---------------------------------------------------------------------------
if nargin < 2
    order = 'T*R*Z*S';
end;

T  =   [1   0   0   P(1);
        0   1   0   P(2);
        0   0   1   P(3);
        0   0   0   1];

R1  =  [1    0      0          0;
        0    cos(P(4))  sin(P(4))  0;
        0   -sin(P(4))  cos(P(4))  0;
        0    0      0          1];

R2  =  [cos(P(5))  0    sin(P(5))  0;
        0          1    0      0;
       -sin(P(5))  0    cos(P(5))  0;
        0          0    0          1];

R3  =  [cos(P(6))   sin(P(6))   0  0;
       -sin(P(6))   cos(P(6))   0  0;
        0           0           1  0;
        0           0       0  1];

R   = R1*R2*R3;

Z   =  [P(7)    0       0       0;
        0       P(8)    0       0;
        0       0       P(9)    0;
        0       0       0       1];

S   =  [1       P(10)   P(11)   0;
        0       1   P(12)   0;
        0       0       1   0;
        0       0       0       1];

A = eval(sprintf('%s;', order));
if ~isnumeric(A) || ndims(A) ~= 2 || any(size(A) ~= 4)
    error('Order expression ''%s'' did not return a valid 4x4 matrix.', ...
          order);
end;