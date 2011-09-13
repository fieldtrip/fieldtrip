function vol = ft_headmodel_bem_openmeeg(geom, varargin)

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
%   vol = ft_headmodel_openmeeg(geom, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   isolatedsource   = string, 'yes' or 'no'
%   hdmfile          = string, filename with BEM headmodel
%   conductivity     = vector, conductivity of each compartment
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

ft_hastoolbox('openmeeg', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the first part is largely shared with the dipoli and bemcp implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the optional arguments
isolatedsource  = keyval('isolatedsource', varargin);
hdmfile         = keyval('hdmfile', varargin);
conductivity    = keyval('conductivity', varargin);

% start with an empty volume conductor
vol = [];

if ~isempty(hdmfile)
  hdm = ft_read_vol(hdmfile);
  % copy the boundary of the head model file into the volume conduction model
  vol.bnd = hdm.bnd;
  if isfield(hdm, 'cond')
    % also copy the conductivities
    vol.cond = hdm.cond;
  end
else
  % copy the boundaries from the geometry into the volume conduction model
  vol.bnd = geom.bnd;
end

% determine the number of compartments
numboundaries = length(vol.bnd);

if isempty(isolatedsource)
  if numboundaries>1
    % the isolated source compartment is by default the most inner one
    isolatedsource = true;
  else
    isolatedsource = false;
  end
else
  % convert into a boolean
  isolatedsource = istrue(cfg.isolatedsource);
end

if ~isfield(vol, 'cond')
  % assign the conductivity of each compartment
  vol.cond = conductivity;
end

% determine the nesting of the compartments
nesting = zeros(numboundaries);
for i=1:numboundaries
  for j=1:numboundaries
    if i~=j
      % determine for a single vertex on each surface if it is inside or outside the other surfaces
      curpos = vol.bnd(i).pnt(1,:); % any point on the boundary is ok
      curpnt = vol.bnd(j).pnt;
      curtri = vol.bnd(j).tri;
      nesting(i,j) = bounding_mesh(curpos, curpnt, curtri);
    end
  end
end

if sum(nesting(:))~=(numboundaries*(numboundaries-1)/2)
  error('the compartment nesting cannot be determined');
end

% for a three compartment model, the nesting matrix should look like
%    0 0 0     the first is the most outside, i.e. the skin
%    0 0 1     the second is nested inside the 3rd, i.e. the outer skull
%    0 1 1     the third is nested inside the 2nd and 3rd, i.e. the inner skull
[~, order] = sort(sum(nesting,2));

fprintf('reordering the boundaries to: ');
fprintf('%d ', order);
fprintf('\n');

% update the order of the compartments
vol.bnd    = vol.bnd(order);
vol.cond   = vol.cond(order);
vol.skin_surface   = 1;
vol.source = numboundaries;

if isolatedsource
  fprintf('using compartment %d for the isolated source approach\n', vol.source);
else
  fprintf('not using the isolated source approach\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this uses an implementation that was contributed by INRIA Odyssee Team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show the license once
openmeeg_license

% check that the binaries are ok
om_checkombin;

if vol.skin_surface ~= 1
  error('the first compartment should be the skin, the last the source');
end

% flip faces for openmeeg convention
for ii=1:length(vol.bnd)
  vol.bnd(ii).tri = fliplr(vol.bnd(ii).tri);
end

% store the current path and change folder to the temporary one
tmpfolder = cd;

try
  cd(tempdir)
  
  % write the triangulations to file
  bndfile = {};
  for ii=1:length(vol.bnd)
    [junk,tname] = fileparts(tempname);
    bndfile{ii} = [tname '.tri'];
    om_save_tri(bndfile{ii}, vol.bnd(ii).pnt, vol.bnd(ii).tri);
  end
  
  % these will hold the shell script and the inverted system matrix
  [~,tname] = fileparts(tempname);
  if ~ispc
    exefile = [tname '.sh'];
  else
    exefile = [tname '.bat'];
  end
  
  [~,tname] = fileparts(tempname);
  condfile  = [tname '.cond'];
  [~,tname] = fileparts(tempname);
  geomfile  = [tname '.geom'];
  [~,tname] = fileparts(tempname);
  hmfile    = [tname '.bin'];
  [~,tname] = fileparts(tempname);
  hminvfile = [tname '.bin'];
  
  % write conductivity and geometry files
  om_write_geom(geomfile,bndfile);
  om_write_cond(condfile,vol.cond);
  
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
  vol.mat = om_load_sym(hminvfile,'binary');
  cleaner(vol,bndfile,condfile,geomfile,hmfile,hminvfile,exefile)
  cd(tmpfolder)
catch
  warning('an error ocurred while running OpenMEEG');
  disp(lasterr);
  cleaner(vol,bndfile,condfile,geomfile,hmfile,hminvfile,exefile)
  cd(tmpfolder)
end

% remember the type of volume conduction model
vol.type = 'openmeeg';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleaner(vol,bndfile,condfile,geomfile,hmfile,hminvfile,exefile)
% delete the temporary files
for i=1:length(vol.bnd)
  delete(bndfile{i})
end
delete(condfile);
delete(geomfile);
delete(hmfile);
delete(hminvfile);
delete(exefile);

