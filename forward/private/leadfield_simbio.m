function [lf] = leadfield_simbio(dip, elc, vol)

% LEADFIELD_SIMBIO leadfields for a set of dipoles
%
% [lf] = leadfield_simbio(dip, elc, vol);
%
% with input arguments
%   dip     positions of the dipoles (matrix of dimensions NX3)
%   elc     positions of the electrodes (matrix of dimensions MX3)
% 
% and vol being a structure with the element
%   vol.meshfile file containing the 3D mesh filename
%
% The output lf is the leadfields matrix of dimensions M (rows) X N*3 (cols)

% Copyright (C) 2011, Cristiano Micheli

% store the current path and change folder to the temporary one
tmpfolder = cd;

try
  if ~ispc
    cd(tempdir)
    [junk,tname] = fileparts(tempname);
    exefile = [tname '.sh'];
    [junk,tname] = fileparts(tempname);
    elcfile  = [tname '.elc'];
    [junk,tname] = fileparts(tempname);
    dipfile  = [tname '.dip'];
    [junk,tname] = fileparts(tempname);
    outfile  = [tname];
    
    % write the electrodes and dipoles positions in the temporary folder
    sb_write_elc(elc.pnt,elc.labels,elcfile);
    sb_write_dip(dip,dipfile);
    
    % Exe file
    % FIXME: does SimBio have a switch for parallel processing (to run on more cores)?
    efid = fopen(exefile, 'w');
    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['ipm_linux_opt_Venant -i sourcesimulation -h ' vol.meshfile ' -s ./' elcfile, ...
        ' -dip ./' dipfile ' -o ' outfile ' -p parameterfile.par -fwd FEM -sens EEG 2>&1 > /dev/null\n']);
    end
    fclose(efid);
    
    dos(sprintf('chmod +x %s', exefile));
    dos(['./' exefile]);
    cleaner(dipfile,elcfile,outfile,exefile)
    [lf] = sb_read_msr(outfile);
    cd(tmpfolder)
  end

catch
  warning('an error occurred while running SimBio');
  rethrow(lasterror)
  cleaner(dipfile,elcfile,outfile,exefile)
  cd(tmpfolder)
end

function cleaner(dipfile,elcfile,outfile,exefile)
delete(dipfile);
delete(elcfile);
delete(outfile);
delete(exefile);
delete([outfile '.fld']);
delete([outfile '.hex']);
delete([outfile '.pot']);
delete([outfile '.pts']);
delete(['*.v.ascii']);
delete(['*.v.potential']);
