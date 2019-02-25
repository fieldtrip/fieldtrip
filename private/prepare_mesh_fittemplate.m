function M = prepare_mesh_fittemplate(headshape,template)

% PREPARE_MESH_FITTEMPLATE computes an affine transformation matrix between 2 point clouds 
%
% This function relies on cpd toolbox found in the external/cpd folder
%
%
% See also FT_PREPARE_MESH

% Copyright (C) 2019, Simon Homoelle
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% ensure that the input is consistent with what this function expects

% add toolbox cpd
ft_hastoolbox('cpd'); %% has to be edited


% prepare control points
axis_limits = determine_border(template,headshape);
res = 10;
[X, Y, Z]  = ndgrid(linspace(axis_limits(1,1),axis_limits(1,2),res), linspace(axis_limits(2,1),axis_limits(2,2),res), linspace(axis_limits(3,1),axis_limits(3,2),res));
ctrl_pts   = [X(:) Y(:) Z(:)];


config.model        = template;
config.scene        = headshape;
config.ctrl_pts     = ctrl_pts;
config.init_param   = zeros(size(ctrl_pts));
config.init_sigma   = 0.5;
config.anneal_rate  = 0.97;
config.outliers     = 0.5;
config.lambda       = 1;
config.beta         = 1;
config.max_iter     = 50;
config.max_em_iter  = 5;
config.tol          = 1e-18;
config.emtol        = 1e-15;
config.motion       = 'tps';
%config.init_param = zeros(25,2);

[param,model] = gmmreg_cpd(config);
ft_plot_mesh(headshape)
ft_plot_mesh(model,'vertexcolor','red')
% the affine transformation can be found in the beginning of param

syms a b c d e f g h k l m n

eqn = template(1:12,:)*[a b c; d e f;g h k] == model(1:12,:) + repmat([l m n],12,1); 
M = solve(eqn, [a b c d e f g h k l m n]);


end
