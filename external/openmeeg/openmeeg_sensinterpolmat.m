function [gain] = openmeeg_sensinterpolmat(elec, vol)

% openmeeg_sensinterpolmat computes the OpenMEEG sensor interpolation matrix
% adopted from openmeeg_DSM

% store the current path and change folder to the temporary one
tmpfolder = cd;
om_checkombin;

bndom = vol.bnd;

try
    cd(tempdir)

    % write the triangulations to file
    bndfile = {};
    for i=1:length(vol.bnd)
        [~,tname] = fileparts(tempname);
        bndfile{i} = [tname '.tri'];
        ok = checknormals(bndom(i));
        if ~ok
          bndom(i).tri = fliplr(bndom(i).tri);
        end
        om_save_tri(bndfile{i}, bndom(i).pos, bndom(i).tri);
    end
    
    % these will hold the shell script and the inverted system matrix
    [~,tname] = fileparts(tempname);
    if ~ispc
      exefile = [tname '.sh'];
    else
      exefile = [tname '.bat'];
    end

    [~,tname] = fileparts(tempname);
    condfile = [tname '.cond'];
    [~,tname] = fileparts(tempname);
    geomfile = [tname '.geom'];
    [~,tname] = fileparts(tempname);
    electrodefile = [tname '.elec'];
    [~,tname] = fileparts(tempname);
    h2mfile = [tname '.bin'];

    % write conductivity and geometry files
    om_write_geom(geomfile,bndfile);
    om_write_cond(condfile,vol.cond);

    % handle electrode file
    elec=ft_convert_units(elec,vol.unit);
    om_save_full(elec.elecpos,electrodefile,'ascii');

    % Exe file
    efid = fopen(exefile, 'w');
    omp_num_threads = feature('numCores');

    str = ' -H2EM';

    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['export OMP_NUM_THREADS=',num2str(omp_num_threads),'\n']);
      fprintf(efid,['om_assemble' str ' ./',geomfile,' ./',condfile,' ./',electrodefile,' ./',h2mfile,' 2>&1 > /dev/null\n']);
    else
      fprintf(efid,['om_assemble' str ' ./',geomfile,' ./',condfile,' ./',electrodefile,' ./',h2mfile,'\n']);
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
    disp(['Assembling OpenMEEG Sensor Interpolation matrix']);
    stopwatch = tic;
    if ispc
        dos([exefile]);
    else
        dos(['./' exefile]);
    end
    gain = om_load_sparse(h2mfile,'binary');
    toc(stopwatch);
    cleaner(vol,bndfile,condfile,geomfile,exefile,electrodefile,h2mfile)
    cd(tmpfolder)
catch
    warning('an error ocurred while running OpenMEEG');
    disp(lasterr);
    cleaner(vol,bndfile,condfile,geomfile,exefile,electrodefile,h2mfile)
    cd(tmpfolder)
end

function cleaner(vol,bndfile,condfile,geomfile,exefile,electrodefile,h2mfile)
% delete the temporary files
for i=1:length(vol.bnd)
    if exist(bndfile{i},'file'),delete(bndfile{i}),end
end
if exist(condfile,'file'),delete(condfile);end
if exist(geomfile,'file'),delete(geomfile);end
if exist(exefile,'file'),delete(exefile);end
if exist(electrodefile,'file'),delete(electrodefile);end
if exist(h2mfile,'file'),delete(h2mfile);end

function ok = checknormals(bnd)
% FIXME: this method is rigorous only for star shaped surfaces
ok = 0;
pnt = bnd.pos;
tri = bnd.tri;
% translate to the center
org = mean(pnt,1);
pnt(:,1) = pnt(:,1) - org(1);
pnt(:,2) = pnt(:,2) - org(2);
pnt(:,3) = pnt(:,3) - org(3);

w = sum(solid_angle(pnt, tri));

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