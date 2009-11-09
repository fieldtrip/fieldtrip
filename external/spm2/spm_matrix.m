function [A] = spm_matrix(P)
% returns an affine transformation matrix
% FORMAT [A] = spm_matrix(P)
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
% A     - affine transformation matrix
%___________________________________________________________________________
%
% spm_matrix returns a matrix defining an orthogonal linear (translation,
% rotation, scaling or affine) transformation given a vector of
% parameters (P).  The transformations are applied in the following order
% (i.e., the opposite to which they are specified):
%
% 1) shear
% 2) scale
% 3) rotation - yaw, roll & pitch
% 4) translation
%
% SPM uses a PRE-multiplication format i.e. Y = A*X where X and Y are 4 x n
% matrices of n coordinates.
%
%__________________________________________________________________________
% @(#)spm_matrix.m	2.1 02/08/14

% pad P with 'null' parameters
%---------------------------------------------------------------------------
q  = [0 0 0 0 0 0 1 1 1 0 0 0];
P  = [P q((length(P) + 1):12)];

T  =   [1 	0 	0 	P(1);
        0 	1 	0 	P(2);
        0 	0 	1 	P(3);
        0 	0 	0 	1];

R1  =  [1    0   	0   	   0;
        0    cos(P(4))  sin(P(4))  0;
        0   -sin(P(4))  cos(P(4))  0;
        0    0    	0   	   1];

R2  =  [cos(P(5))  0   	sin(P(5))  0;
        0    	   1    0  	   0;
       -sin(P(5))  0  	cos(P(5))  0;
        0          0    0   	   1];

R3  =  [cos(P(6))   sin(P(6))   0  0;
       -sin(P(6))   cos(P(6))   0  0;
        0           0           1  0;
        0     	    0    	0  1];

Z   =  [P(7) 	0   	0    	0;
        0    	P(8) 	0    	0;
        0    	0    	P(9) 	0;
        0    	0    	0    	1];

S   =  [1   	P(10)   P(11)   0;
        0   	1 	P(12)   0;
        0   	0   	1	0;
        0    	0    	0    	1];

A = T*R1*R2*R3*Z*S;
