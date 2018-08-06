function varargout = dartel3(varargin)
% DARTEL 3D image registration stuff
%
%_______________________________________________________________________
%
% FORMAT v = dartel3(v,g,f,param)
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
%                   encodes the means and interpolation of temlate
%                   done using logs and softmax function.
%
% This is for performing a single iteration of the Dartel optimisation.
% All flow fields and images are represented by single precision floating
% point values. Images can be scalar fields, in which case the objective
% function is the sum of squares difference.  Alternatively, images can be
% vector fields, in which case the objective function is the sum of squares
% difference between each scalar field + the sum of squares difference
% between one minus the sum of the scalar fields.
%
%_______________________________________________________________________
%
% FORMAT v = dartel3('cgs',A, b, param)
% v     - the solution
% A     - parameterisation of 2nd derivatives
% b     - parameterisation of first derivatives
% param - 6 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 0 Linear elasticity
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3][4] Voxel sizes
%         - [5][6][7] Regularisation parameters
%           - For "membrane energy", the parameters are
%             lambda, unused and id.
%           - For "linear elasticity", the parameters are
%             mu, lambda, and id
%           - For "bending energy", the parameters are
%             lambda, id1 and id2.
%         - [8] Tolerance.  Indicates required degree of accuracy.
%         - [9] Maximum number of iterations.
%
% This is for solving a set of equations using a conjugate gradient
% solver. This method is less efficient than the Full Multigrid.
% v = inv(A+H)*b
% A, b and v are all single precision floating point.
%
%_______________________________________________________________________
%
% FORMAT v = dartel3('fmg',A, b, param)
% v     - the solution n1*n2*n3*3
% A     - parameterisation of 2nd derivatives 
% b     - parameterisation of first derivatives
% param - 6 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 0 Linear elasticity
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3][4] Voxel sizes
%         - [5][6][7] Regularisation parameters
%           - For "membrane energy", the parameters are
%             lambda, unused and id.
%           - For "linear elasticity", the parameters are
%             mu, lambda, and id
%           - For "bending energy", the parameters are
%             lambda, id1 and id2.
%         - [8] Number of Full Multigrid cycles
%         - [9] Number of relaxation iterations per cycle
%
% Solve equations using a Full Multigrid method.  See Press et al
% for more information.
% v = inv(A+H)*b
% A, b and v are all single precision floating point.
%
%_______________________________________________________________________
%
% FORMAT [y,J] = dartel3('Exp', v, param)
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
% FORMAT m = dartel3('vel2mom', v, param)
% v     - velocity (flow) field n1*n2*n3*3.
% param - 4 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 0 Linear elasticity
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3][4] Voxel sizes
%         - [5][6][7] Regularisation parameters
%           - For "membrane energy", the parameters are
%             lambda, unused and id.
%           - For "linear elasticity", the parameters are
%             mu, lambda, and id
%           - For "bending energy", the parameters are
%             lambda, id1 and id2.
% m       - `momentum' field n1*n2*n3*3.
%
% Convert a flow field to a momentum field by m = H*v, where
% H is the large sparse matrix encoding some form of regularisation.
% v and m are single precision floating point.
%
%_______________________________________________________________________
%
% FORMAT y3 = dartel3('comp',y1,y2)
% y1, y2 - deformation fields n1*n2*n3*3.
% y3     - deformation field field n1*n2*n3*3.
%
% Composition of two deformations y3 = y1(y2)
% y1, y2 and y3 are single precision floating point.
%
%
%
% FORMAT [y3,J3] = dartel3('comp', y1, y2, J1, J2)
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
% FORMAT f2 = dartel3('samp', f1, y)
% f1 - input image(s) n1*n2*n3*n4
% y  - points to sample n1*n2*n3*3
% f2 - output image n1*n2*n3*3
%
% Sample a function according to a deformation.
% f2 = f1(y)
% f1, f2 and y are single precision floating point.
%
%_______________________________________________________________________
%
% FORMAT v2 = dartel3('resize', v1, dim)
% v1  - input fields n1*n2*n3*n4
% v2  - output field dim1*dim2*dim3*n4
% dim - output dimensions
%
% Resize a field according to dimensions dim.  This is
% a component of the FMG approach.
%
%_______________________________________________________________________
%
% FORMAT v3 = dartel3('brc', v1, v2)
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
% Note that the boundary conditions are circulant throughout.
% Interpolation is trilinear, except for the resize function
% which uses a 2nd degree B-spline (without first deconvolving).
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: dartel3.m 5506 2013-05-14 17:13:43Z john $

%error('Not compiled for %s in MATLAB %s  (see make.m)\n', computer, version);

if nargin==0,
    error('Incorrect usage.');
elseif ~isa(varargin{1},'char'),
    param0 = varargin{4};
    switch param0(1),
    case 1,
        param = [param0(4) param0(2) 0  0 0];
    case 2,
        param = [param0(4) param0(3) param0(2) 0 0];
    otherwise
        param = [param0(4) 0 0 param0(2) param0(3)];
    end
    param = [param param0(5:end)];
    [varargout{1:nargout}] = spm_diffeo('dartel',varargin{1},varargin{2:3},param);
else
    switch varargin{1},
    case {'fmg','cgs'},
        param0 = varargin{4};
        vx2    = mean(param0(2:4).^2);
        switch param0(1),
        case 1,
            param = [param0(7) 0 0 param0(5) 0  0 0]*vx2;
        case 2,
            param = [param0(7) param0(6) param0(5) 0 0]*vx2;
        otherwise
            param = [param0(7) 0 0 param0(5) param0(6)]*vx2;
        end
        param = [param9(2:4) param param0(8:end)];
        [varargout{1:nargout}] = spm_diffeo('dartel',varargin{1:3},param);

    case {'vel2mom'},
        param0 = varargin{3};
        vx2    = mean(param0(2:4).^2);
        switch param0(1),
        case 1,
            param = [param0(7) param0(5) 0  0 0]*vx2;
        case 2,
            param = [param0(7) param0(6) param0(5) 0 0]*vx2;
        otherwise
            param = [param0(7) 0 0 param0(5) param0(6)]*vx2;
        end
        param = [param0(2:4) param param0(8:end)];
        [varargout{1:nargout}] = spm_diffeo(varargin{1:2},param);

    otherwise
        [varargout{1:nargout}] = spm_diffeo(varargin{:});
    end
end


