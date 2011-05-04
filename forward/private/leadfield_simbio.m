function [lf] = leadfield_simbio(elc, vol, dip)

% LEADFIELD_SIMBIO leadfields for a set of dipoles
%
% [lf] = leadfield_simbio(elc, vol, headmeshfile);
%
% with input arguments
%   elc     positions of the electrodes (matrix of dimensions MX3)
%   vol.brainmesh contains the positions of the vertices on which
%           the potentials are calculated
%   vol.headmesh contains the position of the vertices of the 'head mesh'
%           which is needed to perform a FEM analysis, and an integer
%           number.
%           This is the output of a step in which from a labelled MRI a
%           volumetric mesh of points is generated, and a label (or compartment) 
%           is assigned to each point
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
    parfile  = [tname '.par'];
    [junk,tname] = fileparts(tempname);
    outfile  = [tname];
    
    % write the electrodes and dipoles positions in the temporary folder
    disp('Writing the accessory files on disk...')
    sb_write_elc(elc.pnt,elc.label,elcfile);
    sb_write_dip(dip,dipfile);
    % writes the parameters file
    cfg = [];
    cfg.cond   = vol.cond;
    cfg.labels = vol.labels;
    sb_write_par(cfg,parfile);
    
    % FIXME: vol in future will contain the mesh as well (right?)
    % at the moment contains only a name of the vista file for memory
    % reasons (vol.headmesh)
    
    % Exe file
    % FIXME: does SimBio have a switch for parallel processing (to run on more cores)?
    efid = fopen(exefile, 'w');
    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
%       fprintf(efid,['ipm_linux_opt_Venant -i leadfieldmatrix_onbrainsurface -h ' vol.headmesh ' -s ./' elcfile, ...
%         ' -g ' vol.brainmesh ' -leadfieldfile ' outfile ' -p ' parfile ' -fwd FEM -sens EEG 2>&1 > /dev/null\n']);

      fprintf(efid,['ipm_linux_opt_Venant -i sourcesimulation -h ' vol.headmesh ' -s ./' elcfile, ...
                    ' -dip ' dipfile ' -o ' outfile ' -p ' parfile ' -fwd FEM -sens EEG 2>&1 > /dev/null\n']);
    end
    fclose(efid);
    
    dos(sprintf('chmod +x %s', exefile));
    disp('SimBio is calculating the LeadFields, this may take some time ...')
    
    stopwatch = tic;
    dos(['./' exefile]);
    disp([ 'elapsed time: ' num2str(toc(stopwatch)) ])
    
    [lf] = sb_read_msr(outfile);
    cleaner(dipfile,elcfile,outfile,exefile)
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
