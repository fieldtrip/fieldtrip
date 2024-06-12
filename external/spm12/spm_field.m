function varargout = spm_field(varargin)
% A compiled routine for various spatially regularised inverse problems
%__________________________________________________________________________
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
%__________________________________________________________________________
%
% FORMAT u = spm_field('vel2mom', v, param)
% v     - A field (n1*n2*n3*n4, single).
% param - 6 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6] Regularisation parameters
% u       - Result of applying differential operator (n1*n2*n3*n4, single).
%
% This generates u = A*v, where A is computed as described above.
%__________________________________________________________________________
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
%__________________________________________________________________________
%
% FORMAT spm_field('boundary',b)
% Set the boundary condition.
% b - boundary condition (0 or 1, see above). 
%==========================================================================
%==========================================================================
%
% L1: The following functions are dedicated to L1 types of penalties
%     (total-variation, etc.), when solved using a reweighted least 
%     squares algorithm.
%     Currently, only membrane energy is implemented.
%__________________________________________________________________________
%
% FORMAT u = spm_field('vel2mom1', v, w, param)
% v     - A field (n1*n2*n3*n4, single).
% w     - A field (n1*n2*n3, single) of positive weights.
% param - 4 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4]       Regularisation parameter (membrane energy)
% u       - Result of applying differential operator (n1*n2*n3*n4, single).
%
% This is a generalisation of vel2mom for differential operators that are
% locally weighted. w contains a map of positive weights that are shared
% across channels.
%__________________________________________________________________________
%
% FORMAT u = spm_field('diaginv1', H, w, param)
% H     - Parameterisation of a Hessian at each voxel
%         (n1*n2*n3*(n4*(n4-1)), single)
% w     - A field (n1*n2*n3, single) of positive weights.
% param - 4 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4]       Regularisation parameter (membrane energy)
% u       - diag(inv(H + L)).
%
% This function computes the diagonal of the inverse of the Hessian
% (u = diag(inv(H + L))). To make the inversion tractable, L is 
% approximated by its diagonal. It allows to approximate the posterior
% uncertainty  in a (Bayesian) reweighted least-squares setting. 
%__________________________________________________________________________
%
% FORMAT u = spm_field('trinv1', H, w, param)
% H     - Parameterisation of a Hessian at each voxel
%         (n1*n2*n3*(n4*(n4-1)), single)
% w     - A field (n1*n2*n3, single) of positive weights.
% param - 4 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4]       Regularisation parameter (membrane energy)
% u       - trace(inv(H + L)).
%
% This function computes the trace of the inverse of the Hessian
% (u = trace(inv(H + L))). To make the inversion tractable, L is 
% approximated by its diagonal. It allows to approximate the posterior
% uncertainty  in a (Bayesian) reweighted least-squares setting. 
%__________________________________________________________________________
%
% FORMAT Ap = spm_field('Atimesp', A, p)
% A     - A field of symmetric matrices (n1*n2*n3*(n4*(n4-1)), single)
% p     - A field (n1*n2*n3*n4, single).
% Ap    - A*p.
%
% This function computes efficiently a lot of matrix-vector products.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
error('spm_field.c not compiled - see Makefile')
