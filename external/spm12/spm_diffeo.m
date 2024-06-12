function varargout = spm_diffeo(varargin)
% MEX function called for image registration stuff
%
%__________________________________________________________________________
%
% FORMAT u = spm_diffeo('vel2mom', v, param)
% v     - velocity (flow) field n1*n2*n3*3.
% param - 8 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6][7][8] Regularisation parameters
%           - [4] Absolute displacements need to be penalised by a tiny
%                 amount.  The first element encodes the amount of
%                 penalty on these.  Ideally, absolute displacements
%                 should not be penalised, but it is often necessary
%                 for technical reasons.
%           - [5] The `membrane energy' of the deformation is penalised,
%                 usually by a relatively small amount. This penalises
%                 the sum of squares of the derivatives of the velocity
%                 field (ie the sum of squares of the elements of the
%                 Jacobian tensors).
%           - [6] The `bending energy' is penalised. This penalises the
%                 sum of squares of the 2nd derivatives of the parameters.
%           - [7][8] Linear elasticity regularisation is also included.
%                    The first parameter (mu) is similar to that for
%                    linear elasticity, except it penalises the sum of
%                    squares of the Jacobian tensors after they have been
%                    made symmetric (by averaging with the transpose).
%                    This term essentially penalises length changes,
%                    without penalising rotations.
%                    The final term also relates to linear elasticity,
%                    and is the weight that denotes how much to penalise
%                    changes to the divergence of the velocities (lambda).
%                    This divergence is a measure of the rate of volumetric
%                    expansion or contraction.
% u       - `momentum' field n1*n2*n3*3.
%
% Convert a velocity field to a momentum field by u = A*v, where
% A is the large sparse matrix encoding some form of regularisation.
% v and m are single precision floating point.
%
%__________________________________________________________________________
%
% FORMAT v = spm_diffeo('mom2vel',g, param)
% v     - the solution n1*n2*n3*3
% g     - parameterisation of first derivatives
% param - 10 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6][7][8] Regularisation settings (see vel2mom).
%         - [9] Number of Full Multigrid cycles.
%         - [10] Number of relaxation iterations per cycle.
%
% Solve equations using a Full Multigrid method.  See Press et al
% for more information.
% v = inv(A)*g
% g and v are both single precision floating point.
%
%__________________________________________________________________________
%
% FORMAT v = spm_diffeo('fmg',H, g, param)
% v     - the solution n1*n2*n3*3
% H     - parameterisation of 2nd derivatives 
% g     - parameterisation of first derivatives
% param - 10 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6][7][8] Regularisation settings (see vel2mom).
%         - [9] Number of Full Multigrid cycles.
%         - [10] Number of relaxation iterations per cycle.
%
% Solve equations using a Full Multigrid method, but using Hessian of
% the matching term.  See Press et al for more information.
% v = inv(A+H)*g
% H, g and v are all single precision floating point.
%
%__________________________________________________________________________
%
% FORMAT v = spm_diffeo('cgs',H, g, param)
% v     - the solution
% H     - parameterisation of 2nd derivatives
% g     - parameterisation of first derivatives
% param - 10 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6][7][8] Regularisation settings (see vel2mom).
%         - [9] Tolerance.  Indicates required degree of accuracy.
%         - [10] Maximum number of iterations.
%
% This is for solving a set of equations using a conjugate gradient
% solver. This method is less efficient than the Full Multigrid, and
% is included for illustrative purposes.
% v = inv(A+H)*g
% H, g and v are all single precision floating point.
%
%__________________________________________________________________________
%
% FORMAT F = spm_diffeo('kernel',d,prm)
% d   - image dimensions
% prm - 8 parameters (settings).
%       These are described above (for 'vel2mom').
% F   - The differential operator encoded as an image (or images).
%       Convolving a velocity field by this will give the momentum.
%
%       Note that a smaller (3D) kernel is obtained when the linear
%       elasticity settings are all zero.  If any of the linear
%       elasticity settings are non-zero, the resulting kernel is
%       represented by a 5D array. For the 3D form, the voxel sizes
%       need to be incorporated as an additional scaling of the kernel.
%       See the code in spm_shoot_greens.m for an illustration.
%
%__________________________________________________________________________
%
% FORMAT y3 = spm_diffeo('comp',y1,y2)
% y1, y2 - deformation fields n1*n2*n3*3.
% y3     - deformation field field n1*n2*n3*3.
%
% Composition of two deformations y3 = y1(y2)
% y1, y2 and y3 are single precision floating point.
%
%
% FORMAT [y3,J3] = spm_diffeo('comp', y1, y2, J1, J2)
% y1, y2 - deformation fields n1*n2*n3*3.
% y3     - deformation field n1*n2*n3*3.
% J1, J2 - Jacobian tensor fields n1*n2*n3*3*3.
% J3     - Jacobian tensor field n1*n2*n3*3*3.
%
% Composition of two deformations, with their Jacobian fields.
% All fields are single precision floating point.
%
%__________________________________________________________________________
%
% FORMAT iy = spm_diffeo('invdef',y,d,M1,M2);
%
% iy - inverted deformation field of size d(1)*d(2)*d(3)*3.
% y  - original deformation field.
% M1 - An affine mapping from mm to voxels in the coordinate
%      system of the inverse deformation field.
% M2 - An affine mapping from voxels to mm in the coordinate
%      system of the forward deformation field.
%
% Inversion of a deformation field.
%
% The field is assumed to consist of a piecewise affine transformations,
% whereby each cube jointing 8 neighbouring voxels contains eight
% tetrahedra.  The mapping within each tetrahedron is assumed to be
% affine.
%
%  Reference:
%    J. Ashburner, J. Andersson and K. J. Friston (2000).
%    "Image Registration using a Symmetric Prior - in Three-Dimensions".
%    Human Brain Mapping 9(4):212-225 (appendix).
%__________________________________________________________________________
%
% FORMAT [f,dfx,dfy,dfz] = spm_diffeo('bsplins', c, y,d)
% c          - input image(s) of B-spline coefficients n1*n2*n3*n4
%              - see 'bsplinc'
% y          - points to sample n1*n2*n3*3
% d(1:3)     - degree of B-spline (from 0 to 7) along different dimensions
%              - these must be same as used by 'bsplinc'
% d(4:6)     - 1/0 to indicate wrapping along the dimensions
%
% f           - output image n1*n2*n3*n4
% dfx,dfy,dfz - sampled first derivatives
%
% c, f and y are single precision floating point.
%
% This function takes B-spline basis coefficients from spm_bsplinc,
% and re-convolves them with B-splines centred at the new sample points.
% 
% Note that nearest neighbour interpolation is used instead of 0th
% degree B-splines, and the derivatives of trilinear interpolation are
% returned instead of those of 1st degree B-splines.  The difference is
% extremely subtle.
%
% c, f and y are single precision floating point.
% 
%  References:
%    M. Unser, A. Aldroubi and M. Eden.
%    "B-Spline Signal Processing: Part I-Theory,"
%    IEEE Transactions on Signal Processing 41(2):821-832 (1993).
% 
%    M. Unser, A. Aldroubi and M. Eden.
%    "B-Spline Signal Processing: Part II-Efficient Design and
%    Applications,"
%    IEEE Transactions on Signal Processing 41(2):834-848 (1993).
% 
%    M. Unser.
%    "Splines: A Perfect Fit for Signal and Image Processing,"
%    IEEE Signal Processing Magazine, 16(6):22-38 (1999)
% 
%    P. Thevenaz and T. Blu and M. Unser.
%    "Interpolation Revisited"
%    IEEE Transactions on Medical Imaging 19(7):739-758 (2000).
%
%__________________________________________________________________________
%
% FORMAT c = spm_diffeo('bsplinc',f,d)
%   f - an image
%   d(1:3) - degree of B-spline (from 0 to 7) along different dimensions
%       d(4:6) - 1/0 to indicate wrapping along the dimensions
%   c - returned volume of B-spline coefficients
%
% This function deconvolves B-splines from f, returning
% coefficients, c.  These coefficients are then passed to 'bsplins'
% in order to sample the data using B-spline interpolation.
%
%__________________________________________________________________________
%
% FORMAT f2 = spm_diffeo('samp', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*n4
%
% Sample a function according to a deformation using trilinear interp.
% f2 = f1(y)
% f1, f2 and y are single precision floating point.
% Uses boundary condiditions that wrap around (circulant - identical to
% the 'pullc' option - but retained for backward compatibility).
%
%__________________________________________________________________________
%
% FORMAT f2 = spm_diffeo('pull', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*n4
%
% Sample a function according to a deformation using trilinear interp.
% f2 = f1(y)
% f1, f2 and y are single precision floating point.
% Values sampled outside the field of view of f1 are assigned a value
% of NaN.
%
%__________________________________________________________________________
%
% FORMAT f2 = spm_diffeo('pullc', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*n4
%
% Sample a function according to a deformation using trilinear interp.
% f2 = f1(y)
% f1, f2 and y are single precision floating point.
% Uses boundary condiditions that wrap around (circulant - identical to
% the 'samp' option).
%
%__________________________________________________________________________
%
% FORMAT f2 = spm_diffeo('push', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*n4
%
% Push values of a function according to a deformation.  Note that the
% deformation should be the inverse of the one used with 'samp' or
% 'bsplins'. f1, f2 and y are single precision floating point.
% Voxels in f1 that would be pushed outside the field of view of f2 
% are ignored.
%
%__________________________________________________________________________
%
% FORMAT f2 = spm_diffeo('pushc', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*n4
%
% Push values of a function according to a deformation, but using
% circulant boundary conditions.  Data wraps around (circulant).
% f1, f2 and y are single precision floating point.
%
%__________________________________________________________________________
%
% FORMAT ut = spm_diffeo('pushg', u0, y)
% u0 - input momentum n1*n2*n3*3
% y  - points to sample n1*n2*n3*3
% ut - output momentum n1*n2*n3*3
%
% FORMAT ut = spm_diffeo('pushg', u0, y)
% u0 - input momentum n1*n2*n3*3
% y  - points to sample n1*n2*n3*3
% J  - Jacobian tensor field of y n1*n2*n3*3*3
% ut - output momentum n1*n2*n3*3
%
% Push values of a momentum field according to a deformation using
% circulant boundary conditions.  This essentially computes
% (Ad_y)^* u = |det dy| (dy)^T u(y), which is a key to the
% EPdiff equations used for geodesic shooting.
% u0, ut and y are single precision floating point.
%
%__________________________________________________________________________
%
% FORMAT f2 = spm_diffeo('resize', f1, dim)
% f1  - input fields n1*n2*n3*n4
% f2  - output field dim1*dim2*dim3*n4
% dim - output dimensions
%
% Resize a field according to dimensions dim.  This is a component of
% the multigrid approach, and is used for prolongation.
%
%__________________________________________________________________________
%
% FORMAT v2 = spm_diffeo('restrict', v1)
% v1  - input fields n1*n2*n3*n4
% v2  - output field dim1*dim2*dim3*n4
%
% Restricts a field such that its dimensions are approximately half
% their original.  This is a component of the multigrid approach.
%
%__________________________________________________________________________
%
% FORMAT J = spm_diffeo('def2jac',y)
% y - Deformation field
% J - Jacobian tensor field of y
%
% Compute Jacobian tensors from a deformation.
%
%__________________________________________________________________________
%
% FORMAT J = spm_diffeo('def2det',y)
% y - Deformation field
% j - Jacobian determinant field of y
%
% Compute Jacobian determinants from a deformation.
%
%__________________________________________________________________________
%
% FORMAT j = spm_diffeo('det',J)
% J - Jacobian tensor field
% j - Jacobian determinant field
%
% Compute determinants of Jacobian tensors.
%
%__________________________________________________________________________
%
% FORMAT g = spm_diffeo('grad',v)
% v  - velocity field
% g  - gradient of velocity field
%
% The grad option can be applied to any collection of 3D volumes. If
% the input has dimensions d1 x d2 x d3 x d4 x d5..., then the output
% has dimensions d1 x d2 x d3 x (d4xd5...) x 3.
%
%__________________________________________________________________________
%
% FORMAT dv = spm_diffeo('div',v)
% v  - velocity field
% dv - divergences of velocity field
%
% Computes divergence from velocity field.  This is indicative of rates
% of volumetric expansion/contraction.
%
%__________________________________________________________________________
%
% FORMAT [y,J] = spm_diffeo('smalldef',v,s)
% v - velocity field
% s - scaling factor
% y - small deformation
% J - approximate Jacobian tensors of small deformation (computed via
%     a matrix exponsntial of the Jacobians of the velocity field).
%
% This function is used for each time step of geodesic shooting.  It may
% change in future to use some form of Pade approximation of the
% small deformation.
%
%__________________________________________________________________________
%
% FORMAT t = spm_diffeo('trapprox',H, param)
% v     - the solution n1*n2*n3*3
% H     - parameterisation of 2nd derivatives 
% param - 10 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6][7][8] Regularisation settings (see vel2mom).
% t     - approximation of [trace((L+H)\L) trace((L+H)\H)];
%
% Generate an approximation of Trace((L+H)\L) and Trace((L+H)\H) for
% to give a ball-park figure for the "degrees of freedom" in Laplace
% approximations.  L is the regulariser in sparse matrix form.  The
% approximation is a poor one, which assumes all the off-diagonals of L
% are 0.
% H is single precision floating point.
%
%__________________________________________________________________________
%
% FORMAT v = spm_diffeo('dartel',v,g,f,param)
% v     - flow field n1*n2*n3*3 (single precision float)
% g     - first image n1*n2*n3*n4 (single precision float)
% f     - second image n1*n2*n3*n4 (single precision float)
% param - 9 parameters (settings)
%       - [1][2][3][4][5] Regularisation parameters
%       - [1] Absolute displacements need to be penalised by a tiny
%             amount.  The first element encodes the amount of
%             penalty on these.  Ideally, absolute displacements
%             should not be penalised.
%       - [2] The `membrane energy' of the deformation is penalised,
%             usually by a relatively small amount. This penalises
%             the sum of squares of the derivatives of the velocity
%             field (ie the sum of squares of the elements of the
%             Jacobian tensors).
%       - [3] The `bending energy' is penalised. This penalises the
%             sum of squares of the 2nd derivatives of the velocity.
%       - [4][5] Linear elasticity regularisation is also included.
%                The first parameter (mu) is similar to that for
%                linear elasticity, except it penalises the sum of
%                squares of the Jacobian tensors after they have been
%                made symmetric (by averaging with the transpose).
%                This term essentially penalises length changes,
%                without penalising rotations.
%                The final term also relates to linear elasticity,
%                and is the weight that denotes how much to penalise
%                changes to the divergence of the velocities (lambda).
%                This divergence is a measure of the rate of volumetric
%                expansion or contraction.
%         - [6] Number of Full Multigrid cycles
%         - [7] Number of relaxation iterations per cycle
%         - [8] K, such that 2^K time points are used to
%               generate the deformations.  A value of zero
%               indicates a small deformation model.
%         - [9] code of 0, 1 or 2.
%               0 - asymmetric sums of squares objective function.
%               1 -  symmetric sums of squares objective function.
%               2 - assumes multinomial distribution, where template
%                   encodes the means and interpolation of template
%                   done using logs and softmax function.
%
% This is for performing a single iteration of the Dartel optimisation.
% All velocity fields and images are represented by single precision floating
% point values. Images can be scalar fields, in which case the objective
% function is the sum of squares difference.  Alternatively, images can be
% vector fields, in which case the objective function is the sum of squares
% difference between each scalar field + the sum of squares difference
% between one minus the sum of the scalar fields.
%
%__________________________________________________________________________
%
% FORMAT [y,J] = spm_diffeo('Exp', v, param)
% v - flow field
% J - Jacobian. Usually a tensor field of Jacobian matrices, but can
%     be a field of Jacobian determinants.
% param - 2 (or 3) parameters.
%         [1] K, the number of recursions (squaring steps), such
%             that exponentiation is done using an Euler-like
%             integration with 2^K time steps.
%         [2] a scaling parameter.
%         If there is a third parameter, and it is set to 1, then
%         the J will be the Jacobian determinants.
%
% A flow field is "exponentiated" to generate a deformation field
% using a scaling and squaring approach.  See the work of Arsigny
% et al, or Cleve Moler's "19 Dubious Ways" papers.
%
%__________________________________________________________________________
%
% Note that the boundary conditions are circulant throughout.
% Interpolation is trilinear, except for the resize function
% which uses a 2nd degree B-spline (without first deconvolving).
%
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
error('spm_diffeo.c not compiled - see Makefile')
