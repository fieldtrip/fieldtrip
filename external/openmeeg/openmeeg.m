function [vol] = openmeeg(vol, isolated)

% OPENMEEG computes a symmetric BEM system matrix
%
% Use as
%   [vol] = openmeeg(vol, isolated)
% 
% Attention: the normals of the mesh describing the volume conductor are by
%  FieldTrip convention pointing outwards (with respect to the mesh center), 
%  whereas OpenMEEG binaries expect them to be poiting inwards.

% Copyright (C) 2009, Alexandre Gramfort
% INRIA Odyssee Project Team


% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information


openmeeg_license
om_checkombin;
sprintf('%s','Calculating BEM model...please wait');

skin   = find_outermost_boundary(vol.bnd);
source = find_innermost_boundary(vol.bnd);

% the first compartment should be the skin, the last the source
% flip the order of the compartments if necessary
if skin==length(vol.bnd) && source==1
    % flip the order of the compartments   
    vol.bnd    = fliplr(vol.bnd(:)');
    vol.skin_surface   = 1;
    vol.source = length(vol.bnd);
elseif skin==1 && source==length(vol.bnd)
    vol.skin_surface   = 1;
    vol.source = length(vol.bnd);
else
    error('the first compartment should be the skin, the last the source');
end

% Flip faces for openmeeg convention
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
