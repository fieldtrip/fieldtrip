function [h2mm, s2mm] = openmeeg_megm(pos, vol, sens)

% OPENMEEG_MEGM computes the OpenMEEG H2MM and S2MM matrices, 
% i.e. the contribution to MEG from the potential and from the source
%
% Use as
%   [h2mm,s2mm] = openmeeg_megm(vol, isolated)

% Copyright (C) 2010-2017, OpenMEEG developers

openmeeg_license;              % show the license (only once)
prefix = om_checkombin;        % check the installation of the binaries

% store the current path and change folder to the temporary one
tmpfolder = cd;

try
    cd(tempdir)

    % write the triangulations to file
    bndfile = {};
    for i=1:length(vol.bnd)
        [junk,tname] = fileparts(tempname);
        bndfile{i} = [tname '.tri'];
        om_save_tri(bndfile{i}, vol.bnd(i).pos, vol.bnd(i).tri);
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
    sqdfile = [tname '.bin'];
    [junk,tname] = fileparts(tempname);
    h2mmfile = [tname '.bin'];
    [junk,tname] = fileparts(tempname);
    s2mmfile = [tname '.bin'];

    % write conductivity and geometry files
    om_write_geom(geomfile,bndfile);
    om_write_cond(condfile,vol.cond);

    % handle dipole file
    ndip = size(pos,1);
    pos = [kron(pos,ones(3,1)),kron(ones(ndip,1),eye(3))]; % save pos with each 3D orientation
    om_save_full(pos,dipfile,'ascii');

    % handle squids file
    om_save_full([sens.coilpos,sens.coilori],sqdfile,'ascii');

    % Exe file
    efid = fopen(exefile, 'w');
    omp_num_threads = feature('numCores');

    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['export OMP_NUM_THREADS=',num2str(omp_num_threads),'\n']);
      fprintf(efid,[prefix 'om_assemble -DS2MM ./',dipfile,' ./',sqdfile,' ./',s2mmfile,' 2>&1 > /dev/null\n']);
      fprintf(efid,[prefix 'om_assemble -H2MM ./',geomfile,' ./',condfile,' ./',sqdfile,' ./', h2mmfile,' 2>&1 > /dev/null\n']);
    else
      fprintf(efid,[prefix 'om_assemble -DS2MM ./',dipfile,' ./',sqdfile,' ./',s2mmfile,'\n']);
      fprintf(efid,[prefix 'om_assemble -H2MM ./',geomfile,' ./',condfile,' ./',sqdfile,' ./', h2mmfile,'\n']);
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
    disp(['Assembling OpenMEEG H2MM and S2MM matrices']);
    stopwatch = tic;
    if ispc
        dos(exefile);
    else
        dos(['./' exefile]);
    end
    h2mm = om_load_full(h2mmfile,'binary');
    s2mm = om_load_full(s2mmfile,'binary');
    toc(stopwatch);
    cleaner(vol,bndfile,condfile,geomfile,exefile,dipfile,h2mmfile,s2mmfile,sqdfile)
    cd(tmpfolder)
catch
    warning('an error ocurred while running OpenMEEG');
    disp(lasterr);
    cleaner(vol,bndfile,condfile,geomfile,exefile,dipfile,h2mmfile,s2mmfile,sqdfile)
    cd(tmpfolder)
end

function cleaner(vol,bndfile,condfile,geomfile,exefile,dipfile,h2mmfile,s2mmfile,sqdfile)
  % delete the temporary files
  for i=1:length(vol.bnd)
      if exist(bndfile{i},'file'),delete(bndfile{i}),end
  end
  if exist(condfile,'file'),delete(condfile);end
  if exist(geomfile,'file'),delete(geomfile);end
  if exist(exefile,'file'),delete(exefile);end
  if exist(dipfile,'file'),delete(dipfile);end
  if exist(h2mmfile,'file'),delete(h2mmfile);end
  if exist(s2mmfile,'file'),delete(s2mmfile);end
  if exist(sqdfile,'file'),delete(sqdfile);end

