function headmodel = ft_headmodel_openmeeg(mesh, varargin)

% FT_HEADMODEL_OPENMEEG creates a volume conduction model of the head using the
% boundary element method (BEM). This function takes as input the triangulated
% surfaces that describe the boundaries and returns as output a volume conduction
% model which can be used to compute leadfields.
%
% This function implements
%   Gramfort et al. OpenMEEG: opensource software for quasistatic
%   bioelectromagnetics. Biomedical engineering online (2010) vol. 9 (1) pp. 45
%   http://www.biomedical-engineering-online.com/content/9/1/45
%   doi:10.1186/1475-925X-9-45
% and
%   Kybic et al. Generalized head models for MEG/EEG: boundary element method
%   beyond nested volumes. Phys. Med. Biol. (2006) vol. 51 pp. 1333-1346
%   doi:10.1088/0031-9155/51/5/021
%
% This link with FieldTrip is derived from the OpenMEEG project with contributions
% from Daniel Wong and Sarang Dalal, and uses external command-line executables.
% See http://openmeeg.github.io/
%
% Use as
%   headmodel = ft_headmodel_openmeeg(bnd, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = vector, conductivity of each compartment
%   tissue           = cell-array with the tissue labels for each compartment
%   checkmesh        = 'yes' or 'no'
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2010-2024, Robert Oostenveld
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

ft_hastoolbox('openmeeg', 1);  % add to path (if not yet on path)
openmeeg_license;              % show the license (only once)
prefix = om_checkombin;        % check the installation of the binaries
if(~ispc) % if Linux/Mac, set number of threads
  omp_num_threads = feature('numCores');
  prefix = ['export OMP_NUM_THREADS=' num2str(omp_num_threads) ' && ' prefix];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the first part is largely shared with the dipoli and bemcp implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the optional arguments
conductivity    = ft_getopt(varargin, 'conductivity');
tissue          = ft_getopt(varargin, 'tissue');
checkmesh       = ft_getopt(varargin, 'checkmesh', 'yes');

% convert to Boolean value
checkmesh = istrue(checkmesh);

% copy the boundaries from the mesh into the volume conduction model
if isfield(mesh, 'bnd')
  mesh = mesh.bnd;
end

if ischar(tissue)
  % it should be a cell-array
  tissue = {tissue};
end

% rename pnt into pos
mesh = fixpos(mesh);

% determine the number of compartments
numboundaries = length(mesh);

% OpenMEEG v2.3 and up internally adjusts the convention for surface
% normals, but OpenMEEG v2.2 expects surface normals to point inwards;
% this checks and corrects if needed
for i=1:numboundaries
  switch surface_orientation(mesh(i).pos, mesh(i).tri)
    case 'outward'
      ft_warning('flipping mesh %d', i);
      mesh(i).tri = fliplr(mesh(i).tri);
    case 'inward'
      % this is ok
    case 'otherwise'
      ft_error('incorrect mesh %d', i)
  end
end

% determine the desired nesting of the compartments
order = surface_nesting(mesh, 'outsidefirst');

% reorder the boundaries
if numel(mesh)>1
  fprintf('reordering the boundaries to: ');
  fprintf('%d ', order);
  fprintf('\n');
  % update the order of the compartments
  mesh = mesh(order);
end

% start with an empty volume conductor
headmodel = [];

if isempty(conductivity)
  ft_warning('no conductivity specified, using default values')
  if numboundaries == 1
    % uniform
    conductivity = 1;
  elseif numboundaries == 3
    % skin/skull/brain
    conductivity = [0.33 0.0042 0.33];
  elseif numboundaries == 4
    % skin/skull/csf/brain
    conductivity = [0.33 0.0042 1.7900 0.33];
  else
    ft_error('conductivity values must be specified')
  end
  headmodel.cond = conductivity;
else
  if numel(conductivity)~=numboundaries
    ft_error('each compartment should have a conductivity value');
  end
  % reorder the user-specified conductivities
  headmodel.cond = conductivity(order);
end

% assign default tissue labels if none provided
% the tissue is used for the filename
if isempty(tissue)
  switch(numboundaries)
    case 1
      tissue = {'scalp'};
    case 2
      tissue = {'skull', 'brain'};
    case 3
      tissue = {'scalp','skull','brain'};
    case 4
      tissue = {'scalp','skull','csf','brain'};
    otherwise
      ft_error('tissue types must be specified')
  end
else
  tissue = tissue(order);
end

% do some sanity checks on the meshes
if checkmesh
  for i=1:numboundaries
    ntri = size(mesh(i).tri,1);
    npos = size(mesh(i).pos,1);
    assert(2*(npos-2)==ntri, 'the number of triangles does not match the number of vertices')
  end

  % check for each of the vertices that it falls inside the previous surface
  for i=2:numboundaries
    for j=1:size(mesh(i).pos,1)
      inside = surface_inside(mesh(i).pos(j,:), mesh(i-1).pos, mesh(i-1).tri);
      assert(inside==true, 'vertex %d of surface %d is outside surface %d', j, i, i-1);
    end
  end

  ft_info('the meshes are closed and properly nested')
end % if checkmesh

headmodel.tissue       = tissue;
headmodel.bnd          = mesh;
headmodel.skin_surface = 1;
headmodel.source       = numboundaries;

% just to be sure
clear mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this uses an implementation that was contributed by INRIA Odyssee Team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

workdir = fullfile(tempdir,['ft_om_' datestr(now,'ddmmyyHHMMSSFFF')]);
mkdir(workdir);

try
  % Write the triangulations to file, named after tissue type.
  bndfile = fullfile(workdir,strcat(tissue,'.tri'));
  for ii=1:length(headmodel.bnd)
    om_save_tri(bndfile{ii}, headmodel.bnd(ii).pos, headmodel.bnd(ii).tri);
  end

  condfile  = fullfile(workdir, 'om.cond');
  geomfile  = fullfile(workdir, 'om.geom');
  hmfile    = fullfile(workdir, 'hm.bin');
  hminvfile = fullfile(workdir, 'hminv.bin');

  % write conductivity and mesh files
  bndlabel = {};
  for i=1:length(headmodel.bnd)
    [dum, bndlabel{i}] = fileparts(bndfile{i});
  end

  om_write_geom(geomfile, bndfile, bndlabel);
  om_write_cond(condfile, headmodel.cond, bndlabel);

  om_status = system([prefix 'om_assemble -HM ' geomfile ' ' condfile ' ' hmfile]);
  if(om_status ~= 0) % status = 0 if successful
    ft_error('Aborting OpenMEEG pipeline due to above error.');
  end

  headmodel.mat = inv(om_load_sym(hmfile,'binary'));

  rmdir(workdir,'s'); % remove workdir with intermediate files

catch me
  rmdir(workdir,'s'); % remove workdir with intermediate files
  rethrow(me);
end

% remember the type of volume conduction model
headmodel.type = 'openmeeg';
