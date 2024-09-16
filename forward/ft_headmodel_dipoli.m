function headmodel = ft_headmodel_dipoli(mesh, varargin)

% FT_HEADMODEL_DIPOLI creates a volume conduction model of the head
% using the boundary element method (BEM) for EEG. This function takes
% as input the triangulated surfaces that describe the boundaries and
% returns as output a volume conduction model which can be used to
% compute leadfields.
%
% This implements
%   Oostendorp TF, van Oosterom A. "Source parameter estimation in
%   inhomogeneous volume conductors of arbitrary shape." IEEE Trans
%   Biomed Eng. 1989 Mar;36(3):382-91.
%
% The implementation of this function uses an external command-line
% executable with the name "dipoli" which is provided by Thom Oostendorp.
%
% Use as
%   headmodel = ft_headmodel_dipoli(mesh, ...)
%
% The mesh is given as a boundary or a struct-array of boundaries (surfaces)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   isolatedsource   = string, 'yes' or 'no'
%   conductivity     = vector, conductivity of each compartment
%   tempdir          = string, allows you to specify the path for the tempory files (default is automatic)
%   tempname         = string, allows you to specify the full tempory name including path (default is automatic)
%   checkmesh        = 'yes' or 'no'
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

ft_hastoolbox('dipoli', 1);

% get the optional arguments
isolatedsource  = ft_getopt(varargin, 'isolatedsource');
conductivity    = ft_getopt(varargin, 'conductivity');
tdir            = ft_getopt(varargin, 'tempdir');
tname           = ft_getopt(varargin, 'tempname');
checkmesh       = ft_getopt(varargin, 'checkmesh', 'yes');

% convert to Boolean value
checkmesh = istrue(checkmesh);

if isfield(mesh, 'bnd')
  mesh = mesh.bnd;
end

% replace pnt with pos
mesh = fixpos(mesh);

% start with an empty volume conductor
headmodel = [];

% determine the number of compartments
numboundaries = numel(mesh);

% Dipoli expects surface normals to point inwards;
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

if isempty(isolatedsource)
  if numboundaries>1
    % the isolated source compartment is by default the most inner one
    isolatedsource = true;
  else
    isolatedsource = false;
  end
else
  % convert into a boolean
  isolatedsource = istrue(isolatedsource);
end

if isolatedsource
  fprintf('using isolated source approach\n');
else
  fprintf('not using isolated source approach\n');
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

if isempty(conductivity)
  ft_warning('no conductivity specified, using default values')
  if numboundaries == 1
    % uniform
    conductivity = 1;
  elseif numboundaries == 3
    % skin/skull/brain
    conductivity = [0.3300 0.0042        0.3300];
  elseif numboundaries == 4
    % skin/skull/csf/brain
    conductivity = [0.3300 0.0042 1.7900 0.3300];
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

headmodel.bnd          = mesh;
headmodel.skin_surface = 1;
headmodel.source       = numboundaries; % this is now the last one

% just to be sure
clear mesh

if isolatedsource
  fprintf('using compartment %d for the isolated source approach\n', headmodel.source);
else
  fprintf('not using the isolated source approach\n');
end

% find the location of the dipoli binary
str = which('dipoli.exe');
[p, f, x] = fileparts(str);
dipoli = fullfile(p, f);  % without the .exe extension
% the following determines the executable to use, which is in all cases the 32-bit version
switch lower(computer)
  case {'maci' 'maci64', 'maca64'}
    % apple computer with 32- or 64-bit Intel processor, or Apple ARM processor (M1 or M2)
    dipoli = [dipoli '.universal'];
  case {'glnx86' 'glnxa64'}
    % linux computer
    dipoli = [dipoli '.glnx86'];
  case {'win32', 'win64', 'pcwin64'}
    % windows computer
    dipoli = [dipoli '.exe'];
  otherwise
    ft_error('there is no dipoli executable for your platform');
end
fprintf('using the executable "%s"\n', dipoli);

% determine the prefix for the temporary files
if isempty(tname)
  if isempty(tdir)
    if ismac
      prefix = tempname('/tmp'); % otherwise the filename with the complete path gets too long
    else
      prefix = tempname;
    end
  else
    prefix = tempname(tdir);
  end
else
  prefix = tname;
end

% write the triangulations to file
bndfile = cell(1,numboundaries);
for i=1:numboundaries
  bndfile{i} = sprintf('%s_%d.tri', prefix, i);
  write_tri(bndfile{i}, headmodel.bnd(i).pos, headmodel.bnd(i).tri);
end


if ispc
  % for now the solution for a PC is slightly different than the solution
  % for the other operating systems: the dipoli.exe works fine with -g 
  % arguments, which specify the boundaries. on the other hand, the other
  % OSs do not seem to (anecdotally) swallow this. therefore, for the other
  % operating systems, the 'interactive' mode is mimicked, by means of the <<EOF syntax

  % these will hold the shell script and the inverted system matrix
  exefile = [prefix '.bat'];
  amafile = [prefix '.ama'];

  fid = fopen(exefile, 'w');

  fprintf(fid, '%s -i %s ', dipoli, amafile);
  for i=1:numboundaries
    fprintf(fid, '-g ');
    if isolatedsource && headmodel.source==i
      % the isolated potential approach should be applied using this compartment
      % the dipoli.exe -h mentions a '%' sign, but Thom told me that a '!' also 
      % should work, and anecdotally this works more robustly, depending on the 
      % character with which the path starts
      fprintf(fid, '!');
      fprintf(fid, '%s ', bndfile{i});
    else
      fprintf(fid, '%s ', bndfile{i});
    end
    fprintf(fid, '%g ', headmodel.cond(i));
  end
  fclose(fid);
else
  % these will hold the shell script and the inverted system matrix
  exefile = [prefix '.sh'];
  amafile = [prefix '.ama'];

  fid = fopen(exefile, 'w');
  fprintf(fid, '#!/bin/sh\n');
  fprintf(fid, '\n');
  fprintf(fid, '%s -i %s << EOF\n', dipoli, amafile);
  for i=1:numboundaries
    if isolatedsource && headmodel.source==i
      % the isolated potential approach should be applied using this compartment
      fprintf(fid, '!%s\n', bndfile{i});
    else
      fprintf(fid, '%s\n', bndfile{i});
    end
    fprintf(fid, '%g\n', headmodel.cond(i));
  end
  fprintf(fid, '\n');
  fprintf(fid, '\n');
  fprintf(fid, 'EOF\n');
  fclose(fid);
  % ensure that the temporary shell script can be executed
  dos(sprintf('chmod +x %s', exefile));
end

try
  % execute dipoli and read the resulting file
  dos(exefile);
  ama = loadama(amafile);
  headmodel = ama2headmodel(ama);

catch
  ft_error('an error ocurred while running the dipoli executable - please look at the screen output');
end

% delete the temporary files
for i=1:numboundaries
  delete(bndfile{i})
end
delete(amafile);
delete(exefile);

% maintain the general FieldTrip convention of outward-oriented surfaces
% (also for visualization). The calculation of the head model was done with
% inward-oriented surfaces (sometimes flipped from the original input).
for i=1:numboundaries
  headmodel.bnd(i).tri = fliplr(headmodel.bnd(i).tri);
end

% remember that it is a dipoli model
headmodel.type = 'dipoli';
