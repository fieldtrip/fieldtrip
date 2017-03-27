function varargout = spm_field(varargin)
% A compiled routine for various spatially regularised inverse problems.
%_______________________________________________________________________
%
% FORMAT v = spm_field(H, g, param)
% v     - the solution (n1*n2*n3*n4, single)
% H     - parameterisation of a Hessian at each voxel
%         (n1*n2*n3*(n4*(n4-1)), single)
%         Because the Hessian is symmetric, elements along the
%         4th dimension are ordered:
%         h(1,1), h(2,2), h(3,3),... h(1,2), h(1,3), ..., h(2,3)...
%         Each vector along the 4th dimension should encode a
%         positive (semi)definite matrix.
% g     - parameterisation of first derivatives (n1*n2*n3*n4, single)
% param - 10 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6] Regularisation settings.
%           - [4] Penalty on absolute values.
%           - [5] Penalty on the `membrane energy'. This penalises
%                 the sum of squares of the gradients of the values.
%           - [6] Penalty on the `bending energy'. This penalises
%                 the sum of squares of the 2nd derivatives.
%         - [7]       Number of Full Multigrid cycles.
%         - [8]       Number of relaxation iterations per cycle.
%
% The function solves equations using a Full Multigrid method (see
% Press et al for more information), but incorporating the Hessian
% of some form of likelihood term.
% v = inv(A+B)*g
%     where A = param(4)*I + param(5)*L + param(6)*L'*L
%     and   I = kron(kron(Iz,Iy),Ix)
%           L = kron(kron(Lz,Iy),Ix) + kron(kron(Iz,Ly),Ix) + kron(kron(Iz,Iy),Lx)
%
%           Ix = eye(n1); Iy = eye(n2); Iz = eye(n3)
%           Lx = toeplitz([2 -1 0 ... 0 -1]/param(1)^2) etc
%
% Note that for ill-conditioned A, some regularisation of the solution
% is included.  This means that the solution is not identical to that
% computed using other methods, it is still appropriate for use in
% Gauss-Newton type optimisation schemes.
% _______________________________________________________________________
%
% FORMAT u = spm_field('vel2mom', v, param)
% v     - A field (n1*n2*n3*n4, single).
% param - 6 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6] Regularisation parameters
% u       - Result of applying differential operator (n1*n2*n3*n4, single).
%
% This generates u = A*v, where A is computed as described above.
% _______________________________________________________________________
%
% FORMAT b = spm_field('boundary')
% Get the current boundary condition.
% b - boundary condition
%     0 - field wraps around at the boundary, as if the field is on a
%         torus (circulant).  This is typically assumed when using
%         FFTs for convolution etc.
%     1 - Neumann boundary condition.
%     Note that after a `clear functions' in MATLAB, the boundary
%     condition is reset to 0.
% _______________________________________________________________________
%
% FORMAT spm_field('boundary',b)
% Set the boundary condition.
% b - boundary condition (0 or 1, see above). 
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_field.m 5981 2014-05-13 12:47:14Z john $


%-This is merely the help file for the compiled routine
error('spm_field.c not compiled - see Makefile')
