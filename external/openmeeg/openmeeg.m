function [vol] = openmeeg(vol, isolated)

% OPENMEEG computes a symmetric BEM system matrix
%
% Use as
%   [vol] = openmeeg(vol, isolated)

% Copyright (C) 2009, Alexandre Gramfort
% INRIA Odyssee Project Team

%$Log: openmeeg.m,v $
%Revision 1.11  2009/10/28 14:55:26  crimic
%added wait message
%
%Revision 1.10  2009/09/24 13:40:52  crimic
%added check to gcc version for linux
%
%Revision 1.9  2009/09/22 14:57:26  alegra
%adding multicore option with env variable and fix windows pb when running exefile
%
%Revision 1.8  2009/08/01 11:51:22  alegra
%adding license and bibtex citation
%
%Revision 1.7  2009/07/24 09:48:27  alegra
%fix pb with tmp files and the change of directory when running the shell scripts
%
%Revision 1.6  2009/05/29 12:51:33  crimic
%minor changes
%
%Revision 1.5  2009/05/11 17:30:11  crimic
%upgraded to meet windows and linux compatibility
%

openmeeg_license
sprintf('%s','Calculating BEM model...please wait')

skin   = find_outermost_boundary(vol.bnd);
source = find_innermost_boundary(vol.bnd);

% the first compartment should be the skin, the last the source
if skin==1 && source==length(vol.bnd)
    vol.skin   = 1;
    vol.source = length(vol.bnd);
elseif skin==length(vol.bnd) && source==1
    % flip the order of the compartments
    vol.bnd    = fliplr(vol.bnd(:)');
    vol.skin   = 1;
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

% keyboard

try
    % execute OpenMEEG and read the resulting file
    if ispc
        dos([exefile]);
    else
        version = getgccversion;
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

function version = getgccversion
  % checks for gcc compiler version (works if superior to gcc3)
  tmpdir = pwd;
  cd /tmp
  [junk,tname] = fileparts(tempname);
  txtfile  = [tname '.txt'];
  dos(['gcc -v >& ' txtfile]);
  efid = fopen(txtfile);

  tmp = ''; cnt = 1;
  vec = [];
  while ~isnumeric(tmp)
    tmp = fgetl(efid);
    vec{cnt} = tmp;
    cnt = cnt+1;
  end
  fclose(efid);
  delete(txtfile);
  cd(tmpdir);
  tmp = deblank(vec{cnt-2});
  num = findstr('gcc version ',tmp);
  version = str2num(tmp(num+11:num+12));
  
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
return
