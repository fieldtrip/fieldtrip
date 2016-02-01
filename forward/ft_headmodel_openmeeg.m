function headmodel = ft_headmodel_openmeeg(mesh, varargin)

% FT_HEADMODEL_OPENMEEG creates a volume conduction model of the
% head using the boundary element method (BEM). This function takes
% as input the triangulated surfaces that describe the boundaries and
% returns as output a volume conduction model which can be used to
% compute leadfields.
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
% The implementation in this function is derived from the the OpenMEEG project
%  and uses external command-line executables. See http://gforge.inria.fr/projects/openmeeg
% and http://gforge.inria.fr/frs/?group_id=435.
%
% Use as
%   headmodel = ft_headmodel_openmeeg(mesh, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = vector, conductivity of each compartment
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

%$Id$

ft_hastoolbox('openmeeg', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the first part is largely shared with the dipoli and bemcp implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the optional arguments
conductivity    = ft_getopt(varargin, 'conductivity');

% copy the boundaries from the mesh into the volume conduction model
if isfield(mesh,'bnd')
  mesh = mesh.bnd;
end

% start with an empty volume conductor
headmodel = [];
headmodel.bnd = mesh;

% determine the number of compartments
numboundaries = length(headmodel.bnd);

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
  warning('No conductivity is declared, Assuming standard values\n')
  if numboundaries == 1
    conductivity = 1;
  elseif numboundaries == 3
    % skin/skull/brain
    conductivity = [1 1/80 1] * 0.33;
  else
    error('Conductivity values are required for 2 shells. More than 3 shells not allowed')
  end
  headmodel.cond = conductivity;
else
  if numel(conductivity)~=numboundaries
    error('a conductivity value should be specified for each compartment');
  end
  % update the order of the compartments
  headmodel.cond = conductivity(order);
end

headmodel.skin_surface = 1;
headmodel.source = numboundaries;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this uses an implementation that was contributed by INRIA Odyssee Team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show the license once
% openmeeg_license

% check that the binaries are ok
om_checkombin;

% store the current path and change folder to the temporary one
tmpfolder = cd;
bndom = headmodel.bnd;

try
  cd(tempdir)
  
  % write the triangulations to file
  bndfile = {};
  
  for ii=1:length(bndom)
    % check if vertices' normals are inward oriented
    ok = checknormals(bndom(ii));
    if ~ok
      % Flip faces for openmeeg convention (inwards normals)
      fprintf('flipping normals'' direction\n')
      bndom(ii).tri = fliplr(bndom(ii).tri);
    end
  end
   
  for ii=1:length(headmodel.bnd)
    [junk,tname] = fileparts(tempname);
    bndfile{ii} = [tname '.tri'];
    om_save_tri(bndfile{ii}, bndom(ii).pos, bndom(ii).tri);
  end
  
  % these will hold the shell script and the inverted system matrix
  [tmp,tname] = fileparts(tempname);
  if ~ispc
    exefile = [tname '.sh'];
  else
    exefile = [tname '.bat'];
  end
  
  [tmp,tname] = fileparts(tempname);
  condfile  = [tname '.cond'];
  [tmp,tname] = fileparts(tempname);
  geomfile  = [tname '.mesh'];
  [tmp,tname] = fileparts(tempname);
  hmfile    = [tname '.bin'];
  [tmp,tname] = fileparts(tempname);
  hminvfile = [tname '.bin'];
  
  % write conductivity and mesh files
  om_write_geom(geomfile,bndfile);
  om_write_cond(condfile,headmodel.cond);
  
  % Exe file
  efid = fopen(exefile, 'w');
  omp_num_threads = feature('numCores');
  if ~ispc
    fprintf(efid,'#!/usr/bin/env bash\n');
    fprintf(efid,['export OMP_NUM_THREADS=',num2str(omp_num_threads),'\n']);
    fprintf(efid,['om_assemble -HM ./' geomfile ' ./' condfile ' ./' hmfile ' 2>&1 > /dev/null\n']);
    fprintf(efid,['om_minverser ./' hmfile ' ./' hminvfile ' 2>&1 > /dev/null\n']);
  else
    fprintf(efid,['om_assemble -HM ./' geomfile ' ./' condfile ' ./' hmfile '\n']);
    fprintf(efid,['om_minverser ./' hmfile ' ./' hminvfile '\n']);
  end
  
  fclose(efid);
  
  if ~ispc
    dos(sprintf('chmod +x %s', exefile));
  end
catch
  cd(tmpfolder)
  rethrow(lasterror)
end

try
  % execute OpenMEEG and read the resulting file
  if ispc
    dos([exefile]);
  else
    version = om_getgccversion;
    if version>3
      dos(['./' exefile]);
    else
      error('non suitable GCC compiler version (must be superior to gcc3)');
    end
  end
  headmodel.mat = om_load_sym(hminvfile,'binary');
  cleaner(headmodel,bndfile,condfile,geomfile,hmfile,hminvfile,exefile)
  cd(tmpfolder)
catch
  warning('an error ocurred while running OpenMEEG');
  disp(lasterr);
  cleaner(headmodel,bndfile,condfile,geomfile,hmfile,hminvfile,exefile)
  cd(tmpfolder)
end

% remember the type of volume conduction model
headmodel.type = 'openmeeg';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleaner(headmodel,bndfile,condfile,geomfile,hmfile,hminvfile,exefile)
% delete the temporary files
for i=1:length(headmodel.bnd)
  delete(bndfile{i})
end
delete(condfile);
delete(geomfile);
delete(hmfile);
delete(hminvfile);
delete(exefile);

function ok = checknormals(bnd)
% FIXME: this method is rigorous only for star shaped surfaces
ok = 0;
pos = bnd.pos;
tri = bnd.tri;
% translate to the center
org = mean(pos,1);
pos(:,1) = pos(:,1) - org(1);
pos(:,2) = pos(:,2) - org(2);
pos(:,3) = pos(:,3) - org(3);

w = sum(solid_angle(pos, tri));

if w<0 && (abs(w)-4*pi)<1000*eps
  ok = 0;
  warning('your normals are outwards oriented\n')
elseif w>0 && (abs(w)-4*pi)<1000*eps
  ok = 1;
%   warning('your normals are inwards oriented')
else
  error('your surface probably is irregular\n')
  ok = 0;
end
