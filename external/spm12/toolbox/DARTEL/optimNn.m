function varargout = optimNn(varargin)
% Full multigrid matrix solver stuff (zero gradient at boundaries)
%_______________________________________________________________________
%
% FORMAT v = optimNn('fmg',A, b, param)
% v     - the solution n1*n2*n3*n4
% A     - parameterisation of 2nd derivatives
%         n1*n2*n3*(n4*(n4+1)/2)
%         The first n4 volumes are the diagonal elements, which are
%         followed by the off-diagonals (note that 2nd derivs are%
%         symmetric).  e.g. if n4=3, then the ordering would be
%         (1,1),(2,2),(3,3),(1,2),(1,3),(2,3)
% b     - parameterisation of first derivatives n1*n2*n3*n4
% param - 6 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3][4] Voxel sizes
%         - [5][6][7] Regularisation parameters
%           - For membrane and bending energy, the parameters
%             are lambda, unused and id.
%         - [8] Number of Full Multigrid cycles
%         - [9] Number of relaxation iterations per cycle
%
%           Note that more cycles and iterations may be needed
%           for bending energy than for membrane energy.
%
% Solve equations using a Full Multigrid method.  See Press et al
% for more information.
% v = inv(A+H)*b
% A, b and v are all single precision floating point.
% H is a large sparse matrix encoded by param(1:7).
% The tensor field encoded by A MUST be positive-definite. If it is not,
% then anything could happen (see references about Fisher scoring for
% help on ensuring that second derivatives are positive definite).
%
%_______________________________________________________________________
%
% FORMAT m = optimNn('vel2mom', v, param)
% v     - velocity (flow) field n1*n2*n3*n4.
% param - 4 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3][4] Voxel sizes
%         - [5][6][7] Regularisation parameters
%           - For membrane and bending energy, the parameters
%             are lambda, unusaed and id.
% m       - `momentum' field n1*n2*n3*n4.
%
% Convert a flow field to a momentum field by m = H*v, where
% H is the large sparse matrix encoding some form of regularisation.
% v and m are single precision floating point. This function has uses
% beyond only image registration.
%
%_______________________________________________________________________
%
% Note that the boundary conditions are Neumann (zero gradients at the
% boundaries) throughout. For circulant boundary conditions, use
% use optimN.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: optimNn.m 4758 2012-05-29 15:34:08Z john $

[varargout{1:nargout}] = optim_compat(1,varargin{:});

%error('Not compiled for %s in MATLAB %s  (see make.m)\n', computer, version);

