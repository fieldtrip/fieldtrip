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
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

ft_hastoolbox('dipoli', 1);

% get the optional arguments
isolatedsource  = ft_getopt(varargin, 'isolatedsource');
conductivity    = ft_getopt(varargin, 'conductivity');
tdir            = ft_getopt(varargin, 'tempdir');
tname           = ft_getopt(varargin, 'tempname');

if isfield(mesh, 'bnd')
  mesh = mesh.bnd;
end

% replace pnt with pos
mesh = fixpos(mesh);

% start with an empty volume conductor
headmodel = [];
headmodel.bnd = mesh;

% determine the number of compartments
numboundaries = numel(headmodel.bnd);

% % The following checks can in principle be performed, but are too
% % time-consuming. Instead the code here relies on the calling function to
% % feed in the correct geometry.
% %
% % if ~all(surface_closed(headmodel.bnd))
% %   ft_error('...');
% % end
% % if any(surface_intersection(headmodel.bnd))
% %   ft_error('...');
% % end
% % if any(surface_selfintersection(headmodel.bnd))
% %   ft_error('...');
% % end
%
% % The following checks should always be done.
% headmodel.bnd = surface_orientation(headmodel.bnd, 'outwards'); % might have to be inwards
%
% order = surface_nesting(headmodel.bnd, 'outsidefirst'); % might have to be insidefirst
% headmodel.bnd = headmodel.bnd(order);
%
% FIXME also the cond should be in the right order
%

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
order = surface_nesting(headmodel.bnd, 'outsidefirst');

% rearrange boundaries and conductivities
if numel(headmodel.bnd)>1
  fprintf('reordering the boundaries to: ');
  fprintf('%d ', order);
  fprintf('\n');
  % update the order of the compartments
  headmodel.bnd = headmodel.bnd(order);
end

if isempty(conductivity)
  ft_warning('No conductivity is declared, Assuming standard values\n')
  if numboundaries == 1
    conductivity = 1;
  elseif numboundaries == 3
    % skin/skull/brain
    conductivity = [1 1/80 1] * 0.33;
  elseif numboundaries == 4
    % FIXME: check for better default values here
    % skin / outer skull / inner skull / brain
    conductivity = [1 1/80 1 1] * 0.33;
  else
    ft_error('Conductivity values are required!')
  end
  headmodel.cond = conductivity;
else
  if numel(conductivity)~=numboundaries
    ft_error('a conductivity value should be specified for each compartment');
  end
  headmodel.cond = conductivity(order);
end

headmodel.skin_surface = 1;
headmodel.source       = numboundaries; % this is now the last one

if isolatedsource
  fprintf('using compartment %d for the isolated source approach\n', headmodel.source);
else
  fprintf('not using the isolated source approach\n');
end

% find the location of the dipoli binary
str = which('dipoli.maci');
[p, f, x] = fileparts(str);
dipoli = fullfile(p, f);  % without the .maci extension
% the following determines the executable to use, which is in all cases the 32-bit version
switch lower(computer)
  case {'maci' 'maci64'}
    % apple computer
    dipoli = [dipoli '.maci'];
  case {'glnx86' 'glnxa64'}
    % linux computer
    dipoli = [dipoli '.glnx86'];
  case {'win32', 'win64'}
    % windows computer
    dipoli = [dipoli '.exe'];
  otherwise
    ft_error('there is no dipoli executable for your platform');
end
fprintf('using the executable "%s"\n', dipoli);

% determine the prefix for the temporary files
if isempty(tname)
  if isempty(tdir)
    prefix = tempname;
  else
    prefix = tempname(tdir);
  end
else
  prefix = tname;
end

% checks if normals are inwards oriented, otherwise flips them
for i=1:numel(headmodel.bnd)
  isinward = checknormals(headmodel.bnd(i));
  if ~isinward
    fprintf('flipping the normals inward prior to head matrix calculation\n')
    headmodel.bnd(i).tri = fliplr(headmodel.bnd(i).tri);
  end
end

% write the triangulations to file
bndfile = cell(1,numboundaries);
for i=1:numboundaries
  bndfile{i} = sprintf('%s_%d.tri', prefix, i);
  write_tri(bndfile{i}, headmodel.bnd(i).pos, headmodel.bnd(i).tri);
end

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

try
  % execute dipoli and read the resulting file
  dos(exefile);
  ama = loadama(amafile);
  headmodel = ama2headmodel(ama);
  
catch
  ft_error('an error ocurred while running the dipoli executable - please look at the screen output');
end

% This is to maintain the headmodel.bnd convention (outward oriented), whereas
% in terms of further calculation it should not really matter.
% The calculation fo the head model is done with inward normals
% (sometimes flipped from the original input). This assures that the
% outward oriented mesh is saved outward oriented in the headmodel structure
for i=1:numel(headmodel.bnd)
  isinward = checknormals(headmodel.bnd(i));
  if isinward
    fprintf('flipping the normals outwards after head matrix calculation\n')
    headmodel.bnd(i).tri = fliplr(headmodel.bnd(i).tri);
  end
end

% delete the temporary files
for i=1:numboundaries
  delete(bndfile{i})
end
delete(amafile);
delete(exefile);

% remember that it is a dipoli model
headmodel.type = 'dipoli';


function ok = checknormals(bnd)
% checks if the normals are inward oriented
pos = bnd.pos;
tri = bnd.tri;
% translate to the center
center = median(pos,1);
pos(:,1) = pos(:,1) - center(1);
pos(:,2) = pos(:,2) - center(2);
pos(:,3) = pos(:,3) - center(3);

w = sum(solid_angle(pos, tri));

% FIXME: this method is rigorous only for star shaped surfaces with the origin inside
if w<0 && (abs(w)-4*pi)<1000*eps
  % normals are outwards oriented
  ok = 0;
elseif w>0 && (abs(w)-4*pi)<1000*eps
  % normals are inwards oriented
  ok = 1;
else
  ft_warning('your headmodel boundary is not star-shaped when seen from the center')
  ok = 1;
end
