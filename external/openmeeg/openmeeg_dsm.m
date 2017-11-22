function [dsm] = openmeeg_dsm(pos, vol, flag)

% OPENMEEG_DSM computes the OpenMEEG DSM matrix
%              i.e. Right hand side in the potential equation
%
% Use as
%   [dsm] = openmeeg_dsm(po, vol, flag)
%
% flag = 1 non adaptive algorithm: does not try to approximate the
% potential in the neighborhodd of the leads, by locally refining the BEM surface

% Copyright (C) 2010-2017, OpenMEEG developers

openmeeg_license;              % show the license (only once)
prefix = om_checkombin;        % check the installation of the binaries

% store the current path and change folder to the temporary one
tmpfolder = cd;

bndom = vol.bnd;

try
    cd(tempdir)

    % write the triangulations to file
    bndfile = {};
    for i=1:length(vol.bnd)
        [junk,tname] = fileparts(tempname);
        bndfile{i} = [tname '.tri'];
        ok = checknormals(bndom(i));
        if ~ok
          bndom(i).tri = fliplr(bndom(i).tri);
        end
        om_save_tri(bndfile{i}, bndom(i).pos, bndom(i).tri);
    end
    
    % these will hold the shell script and the inverted system matrix
    [junk,tname] = fileparts(tempname);
    if ~ispc
      exefile = [tname '.sh'];
    else
      exefile = [tname '.bat'];
    end

    [junk,tname] = fileparts(tempname);
    condfile = [tname '.cond'];
    [junk,tname] = fileparts(tempname);
    geomfile = [tname '.geom'];
    [junk,tname] = fileparts(tempname);
    dipfile = [tname '.dip'];
    [junk,tname] = fileparts(tempname);
    dsmfile = [tname '.bin'];

    % write conductivity and geometry files
    om_write_geom(geomfile,bndfile);
    om_write_cond(condfile,vol.cond);

    % handle dipole file
    ndip = size(pos,1);
    pos = [kron(pos,ones(3,1)) , kron(ones(ndip,1),eye(3))]; % save pos with each 3D orientation
    om_save_full(pos,dipfile,'ascii');

    % Exe file
    efid = fopen(exefile, 'w');
    omp_num_threads = feature('numCores');

    if flag
      str = ' -DSMNA';
    else
      str = ' -DSM';
    end
    
    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['export OMP_NUM_THREADS=',num2str(omp_num_threads),'\n']);
      % the following implements Galerkin method and switch can be -DSM or -DSMNA
      % (non adaptive), see OMtrunk/src/assembleSourceMat.cpp, operators.cpp
      fprintf(efid,[prefix 'om_assemble' str ' ./',geomfile,' ./',condfile,' ./',dipfile,' ./',dsmfile,' 2>&1 > /dev/null\n']);
    else
      fprintf(efid,[prefix 'om_assemble' str ' ./',geomfile,' ./',condfile,' ./',dipfile,' ./',dsmfile,'\n']);
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
    disp('Assembling OpenMEEG DSM matrix');
    stopwatch = tic;
    if ispc
        dos(exefile);
    else
        dos(['./' exefile]);
    end
    dsm = om_load_full(dsmfile,'binary');
    toc(stopwatch);
    cleaner(vol,bndfile,condfile,geomfile,exefile,dipfile,dsmfile)
    cd(tmpfolder)
catch
    warning('an error ocurred while running OpenMEEG');
    disp(lasterr);
    cleaner(vol,bndfile,condfile,geomfile,exefile,dipfile,dsmfile)
    cd(tmpfolder)
end

function cleaner(vol,bndfile,condfile,geomfile,exefile,dipfile,dsmfile)
% delete the temporary files
for i=1:length(vol.bnd)
    if exist(bndfile{i},'file'),delete(bndfile{i}),end
end
if exist(condfile,'file'),delete(condfile);end
if exist(geomfile,'file'),delete(geomfile);end
if exist(exefile,'file'),delete(exefile);end
if exist(dipfile,'file'),delete(dipfile);end
if exist(dsmfile,'file'),delete(dsmfile);end

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
%   warning('your normals are outwards oriented\n')
elseif w>0 && (abs(w)-4*pi)<1000*eps
  ok = 1;
%   warning('your normals are inwards oriented')
else
  error('your surface probably is irregular\n')
  ok = 0;
end
