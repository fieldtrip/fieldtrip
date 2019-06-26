function varargout = spm_diffeo(varargin)
% Mex function called for image registration stuff
%
%_______________________________________________________________________
%
% FORMAT u = spm_diffeo('vel2mom', v, param)
% v     - velocity (flow) field n1*n2*n3*3.
% param - 8 parameters (settings)
%         - [1][2][3] Voxel sizes
%         - [4][5][6][7][8] Regularisation parameters
%           - [4] Absolute displacements need to be penalised by a tiny
%                 amount.  The first element encodes the amount of
%                 penalty on these.  Ideally, absolute displacements
%                 should not be penalised, but it is usually necessary
%                 for technical reasons.
%           - [5] The `membrane energy' of the deformation is penalised,
%                 usually by a relatively small amount. This penalises
%                 the sum of squares of the derivatives of the velocity
%                 field (ie the sum of squares of the elements of the
%                 Jacobian tensors).
%           - [6] The `bending energy' is penalised (3rd element). This
%                 penalises the sum of squares of the 2nd derivatives of
%                 the velocity.
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
%_______________________________________________________________________
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
%_______________________________________________________________________
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
%_______________________________________________________________________
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
%_______________________________________________________________________
%
% FORMAT F = spm_diffeo('kernel',d,prm)
% d   - image dimensions
% prm - 8 parameters (settings).
%       These are described above (for 'vel2mom').
% F   - The differential operator encoded as an image (or images).
%       Convolving a velocity field by this will give the momentum.
%
%_______________________________________________________________________
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
%_______________________________________________________________________
%
% FORMAT iy = spm_diffeo('invdef',y,d,M1,M2);
%
% iy - inverted deformation field of size d(1)*d(2)*d(3)*3.
% y  - original deformation field.
% M1 - An affine mapping from mm to voxels in the co-ordinate
%      system of the inverse deformation field.
% M2 - An affine mapping from voxels to mm in the co-ordinate
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
%_______________________________________________________________________
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
% returned insted of those of 1st degree B-splines.  The difference is
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
%_______________________________________________________________________
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
%_______________________________________________________________________
%
% FORMAT f2 = spm_diffeo('samp', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*n4
%
% Sample a function according to a deformation using trilinear interp.
% f2 = f1(y)
% f1, f2 and y are single precision floating point.
%
%_______________________________________________________________________
%
% FORMAT f2 = spm_diffeo('push', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*n4
%
% Push values of a function according to a deformation.  Note that the
% deformation should be the inverse of the one used with 'samp' or 'bsplins'.
% f1, f2 and y are single precision floating point.
%
%_______________________________________________________________________
%
% FORMAT f2 = spm_diffeo('pushc', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*n4
%
% Push values of a function according to a deformation, but using
% circulant boundary conditions.  Data wraps around.
% f1, f2 and y are single precision floating point.
%
%_______________________________________________________________________
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
%_______________________________________________________________________
%
% FORMAT f2 = spm_diffeo('resize', f1, dim)
% f1  - input fields n1*n2*n3*n4
% f2  - output field dim1*dim2*dim3*n4
% dim - output dimensions
%
% Resize a field according to dimensions dim.  This is a component of
% the multigrid approach, and is used for prolongation.
%
%_______________________________________________________________________
%
% FORMAT v2 = spm_diffeo('restrict', v1)
% v1  - input fields n1*n2*n3*n4
% v2  - output field dim1*dim2*dim3*n4
%
% Restricts a field such that its dimensions are approximately half
% their original.  This is a component of the multigrid approach.
%
%_______________________________________________________________________
%
% FORMAT J = spm_diffeo('def2jac',y)
% y - Deformation field
% J - Jacobian tensor field of y
%
% Compute Jacobian tensors from a deformation.
%
%_______________________________________________________________________
%
% FORMAT J = spm_diffeo('def2det',y)
% y - Deformation field
% j - Jacobian determinant field of y
%
% Compute Jacobian determinants from a deformation.
%
%_______________________________________________________________________
%
% FORMAT j = spm_diffeo('det',J)
% J - Jacobian tensor field
% j - Jacobian determinant field
%
% Compute determinants of Jacobian tensors.
%
%_______________________________________________________________________
%
% FORMAT dv = spm_diffeo('div',v)
% v  - velocity field
% dv - divergences of velocity field
%
% Computes divergence from velocity field.  This is indicative of rates
% of volumetric expansion/contraction.
%
%_______________________________________________________________________
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
%_______________________________________________________________________
%
% FORMAT v3 = spm_diffeo('brc', v1, v2)
% v1, v2, v3 - flow fields n1*n2*n3*3
%
% Lie Bracket.  Useful for many things
% e.g. Baker-Campbell-Haussdorf series expansion.
% The Lie bracket is denoted by
% v3 = [v1,v2]
% and on scalar fields, is computed by
% v3 = J1*v2 - J2*v1, where J1 and J2 are the Jacobian
% tensor fields. For matrices, the Lie bracket is simply
% [A,B] = A*B-B*A
%
%_______________________________________________________________________
%
% FORMAT v = spm_diffeo('dartel',v,g,f,param)
% v     - flow field n1*n2*n3*3 (single precision float)
% g     - first image n1*n2*n3*n4 (single precision float)
% f     - second image n1*n2*n3*n4 (single precision float)
% param - 9 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 0 Linear elasticity
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3][4] Regularisation parameters
%           - For "membrane energy", the parameters are
%             lambda, unused and id.
%           - For "linear elasticity", the parameters are
%             mu, lambda, and id
%           - For "bending energy", the parameters are
%             lambda, id1 and id2, such that regularisation is by
%             (-lambda*\grad^2 + id1)^2 + id2
%         - [5] Levenberg-Marquardt regularisation
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
%_______________________________________________________________________
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
%_______________________________________________________________________
%
% Note that the boundary conditions are circulant throughout.
% Interpolation is trilinear, except for the resize function
% which uses a 2nd degree B-spline (without first deconvolving).
%
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_diffeo.m 4890 2012-09-03 15:19:46Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_diffeo.c not compiled - see Makefile')
