function fittemplate = prepare_mesh_fittemplate(cfg,template)

% PREPARE_MESH_FITTEMPLATE creates a individualized template on the basis
% of surface information
%
% This function relies on cpd toolbox found in the external/cpd folder
%
% Configuration options:
%   cfg.method       = 'fittemplate'
%   cfg.headshape    =
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

%
headshape = cfg.headshape.pos;
%% determine outer most layer
index = find_outermost_boundary(template.bnd);
top_template = template.bnd(index).pos;
%% Fit top part

% affine register
opt.corresp = 0;
opt.method  = 'affine';
opt.normalize = 1;
opt.max_it = 100;
opt.fgt=1;
opt.tol = 10e-12;
opt.outliers=0.0;
[transform,~] = cpd_register(headshape,top_template, opt);

% create 4x4 transformation Matrix
M(1:3,1:3)                       = transform.R;
M(1:3,4)                         = transform.t;

%removing structures that became obsolete with changing geometry
fittemplate = template;
if isfield(template,'mat')
    fittemplate = rmfield(fittemplate,'mat');
end    
if isfield(template,'type')
    fittemplate = rmfield(fittemplate,'type');
end    

% warping
for i = 1:length(template.bnd)
fittemplate.bnd(i).pos = ft_warp_apply(M,fittemplate.bnd(i).pos);
end

