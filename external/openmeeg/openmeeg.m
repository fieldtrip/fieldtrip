function [vol] = openmeeg(vol, isolated)

% OPENMEEG computes a symmetric BEM system matrix
%
% Use as
%   [vol] = openmeeg(vol, isolated)
%
% Attention: the normals of the mesh describing the volume conductor are by
%  FieldTrip convention pointing outwards (with respect to the mesh center),
%  whereas OpenMEEG binaries expect them to be poiting inwards.

% Copyright (C) 2010-2017, OpenMEEG developers

warning('OPENMEEG is deprecated, please use FT_PREPARE_HEADMODEL with cfg.method = ''openmeeg'' instead.')

openmeeg_license;              % show the license (only once)
prefix = om_checkombin;        % check the installation of the binaries

sprintf('%s','Calculating BEM model...please wait');

skin   = find_outermost_boundary(vol.bnd);
source = find_innermost_boundary(vol.bnd);

% the first compartment should be the skin, the last the source
% flip the order of the compartments if necessary
if skin==length(vol.bnd) && source==1
  % flip the order of the compartments
  vol.bnd    = fliplr(vol.bnd(:)');
  vol.cond   = fliplr(vol.cond(:)');
  vol.skin_surface   = 1;
  vol.source = length(vol.bnd);
elseif skin==1 && source==length(vol.bnd)
  vol.skin_surface   = 1;
  vol.source = length(vol.bnd);
else
  error('the first compartment should be the skin, the last the source');
end

% store the current path and change folder to the temporary one
tmpfolder = cd;

try
  cd(tempdir)
  % initialize OM boundaries
  bndfile = {};
  bndom = vol.bnd;
  
  % check if normals are outward oriented (as they should be)
  ok = checknormals(bndom);
  
  % Flip faces for openmeeg convention (inwards normals)
  if ~ok
    for ii=1:length(bndom)
      bndom(ii).tri = fliplr(bndom(ii).tri);
    end
  end
  
  % write triangulation files on disk
  for ii=1:length(vol.bnd)
    [junk,tname] = fileparts(tempname);
    bndfile{ii} = [tname '.tri'];
    om_save_tri(bndfile{ii}, bndom(ii).pos, bndom(ii).tri);
  end
  
  % these will hold the shell script and the inverted system matrix
  [junk,tname] = fileparts(tempname);
  if ~ispc
    exefile = [tname '.sh'];
  else
    exefile = [tname '.bat'];
  end
  
  [junk,tname] = fileparts(tempname);
  condfile  = [tname '.cond'];
  [junk,tname] = fileparts(tempname);
  geomfile  = [tname '.geom'];
  [junk,tname] = fileparts(tempname);
  hmfile    = [tname '.bin'];
  [junk,tname] = fileparts(tempname);
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
    fprintf(efid,[prefix 'om_assemble -HM ./' geomfile ' ./' condfile ' ./' hmfile ' 2>&1 > /dev/null\n']);
    fprintf(efid,[prefix 'om_minverser ./' hmfile ' ./' hminvfile ' 2>&1 > /dev/null\n']);
  else
    fprintf(efid,[prefix 'om_assemble -HM ./' geomfile ' ./' condfile ' ./' hmfile '\n']);
    fprintf(efid,[prefix 'om_minverser ./' hmfile ' ./' hminvfile '\n']);
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
    dos(exefile);
  elseif ismac
    dos(['./' exefile]);
  else % assumes linux by default
    version = om_getgccversion;
    if version > 3
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

function cleaner(vol,bndfile,condfile,geomfile,hmfile,hminvfile,exefile)

% delete the temporary files
for ii=1:length(vol.bnd)
  delete(bndfile{ii})
end

delete(condfile);
delete(geomfile);
delete(hmfile);
delete(hminvfile);
delete(exefile);
return

function ok = checknormals(bnd)
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
  % FIXME: this method is rigorous only for star shaped surfaces
  warning('your normals are not oriented correctly')
  ok = 0;
elseif w>0 && abs(w-4*pi)<1000*eps
  ok = 1;
else
  error('your surface probably is irregular')
  ok = 0;
end
