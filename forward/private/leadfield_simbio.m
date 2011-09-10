function [lf] = leadfield_simbio(elc, vol, dip)

% LEADFIELD_SIMBIO leadfields for a set of dipoles
%
% [lf] = leadfield_simbio(elc, vol, headmeshfile);
%
% with input arguments
%   elc     positions of the electrodes (matrix of dimensions MX3)
%           contains the positions of the vertices on which
%           the potentials are calculated
%           There can be 'deep electrodes' too!
%   vol.wf  contains a wireframe structure (or FEM grid) of 'elements'
%           which is needed to perform a FEM analysis.
%   dip     a matrix of dipoles (of dimensions NX3)
%
% The output lf is the leadfields matrix of dimensions M (rows) X N*3 (cols)

% Copyright (C) 2011, Cristiano Micheli

% store the current path and change folder to the temporary one
tmpfolder = cd;

try
  if ~ispc
    cd(tempdir)
    [~,tname] = fileparts(tempname);
    exefile = [tname '.sh'];
    [~,tname] = fileparts(tempname);
    elcfile  = [tname '.elc'];
    [~,tname] = fileparts(tempname);
    dipfile  = [tname '.dip'];    
    [~,tname] = fileparts(tempname);
    parfile  = [tname '.par'];
    [~,tname] = fileparts(tempname);
    outfile  = [tname];
    
    % write the electrodes and dipoles positions in the temporary folder
    disp('Writing the accessory files on disk...')
    
    % FIXME: distinguish between surface and deep electrodes 
    sb_write_elc(elc.pnt,elc.label,elcfile);
    sb_write_dip(dip,dipfile);
    
    % write the parameters file, contains tissue/conductivities match
    % and Simbio call details, mixed together
    sb_write_par(parfile,'cond',vol.cond,'labels',vol.labels);

    % write the vol.wf in a Vista format .v file
    [~,tname] = fileparts(tempname);
    wffile = [tname '.v'];
    % write a vista wireframe file
    write_vista_mesh(wffile,vol.wf.nd,vol.wf.el,vol.labels);
  
    % Exe file
    % FIXME: does SimBio have a switch for parallel processing (to run on more cores)?
    efid = fopen(exefile, 'w');
    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
%       fprintf(efid,['ipm_linux_opt_Venant -i leadfieldmatrix_onbrainsurface -h ' vol.headmesh ' -s ./' elcfile, ...
%         ' -g ' vol.brainmesh ' -leadfieldfile ' outfile ' -p ' parfile ' -fwd FEM -sens EEG 2>&1 > /dev/null\n']);

      fprintf(efid,['ipm_linux_opt_Venant -i sourcesimulation -h ' wffile ' -s ./' elcfile, ...
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
