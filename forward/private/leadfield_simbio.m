function [lf] = leadfield_simbio(dip, elc, vol)

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
  
  % add the Simbio folder to the path
  if ft_hastoolbox('simbio',1,0)
    if ~ispc
      
      cd(tempdir)
      [~,tname] = fileparts(tempname);
      exefile   = [tname '.sh'];
      [~,tname] = fileparts(tempname);
      elcfile   = [tname '.elc'];
      [~,tname] = fileparts(tempname);
      dipfile   = [tname '.dip'];
      [~,tname] = fileparts(tempname);
      parfile   = [tname '.par'];
      [~,tname] = fileparts(tempname);
      meshfile  = [tname '.v'];
      [~,tname] = fileparts(tempname);
      transfermatrix = [tname];
      [~,tname] = fileparts(tempname);
      outfile   = [tname];      
      
      % write the electrodes and dipoles positions in the temporary folder
      disp('Writing the electrodes file on disk...')
      if ~isfield(vol,'deepelec')
        sb_write_elc(warp_apply(inv(vol.transform),elc.pnt),elc.label,elcfile);
      else
        sb_write_elc(warp_apply(inv(vol.transform),elc.pnt),elc.label,elcfile,1);
      end
      disp('Writing the dipoles file on disk...')
      sb_write_dip(warp_apply(inv(vol.transform),dip),dipfile);
      
      % write the parameters file, contains tissues conductivities, the FE
      % grid and Simbio call details, mixed together
      disp('Writing the parameters file on disk...')
      sb_write_par(parfile,'cond',vol.cond,'labels',unique(vol.wf.labels));
      
      % write the vol.wf in a Vista format .v file
      ft_write_headshape(meshfile,vol.wf,'format','vista');
      
      % Exe file
      efid = fopen(exefile, 'w');
      if ~ispc
        fprintf(efid,'#!/usr/bin/env bash\n');
%         fprintf(efid,['ipm_linux_opt_Venant -i FEtransfermatrix -h ./' meshfile ' -s ./' elcfile, ...
%           ' -o ./' transfermatrix ' -p ./' parfile ' -sens EEG 2>&1 > /dev/null\n']); 
%         fprintf(efid,['ipm_linux_opt_Venant -i sourcesimulation -s ./' elcfile ' -t ./' transfermatrix, ...
%           ' -dip ./' dipfile ' -o ./' outfile ' -p ./' parfile ' -fwd FEM -sens EEG 2>&1 > /dev/null\n']);         
        fprintf(efid,['ipm_linux_opt_Venant -i sourcesimulation -h ./' meshfile ' -s ./' elcfile, ...
          ' -dip ./' dipfile ' -o ./' outfile ' -p ./' parfile ' -fwd FEM -sens EEG 2>&1 > /dev/null\n']);
      end
      fclose(efid);
      
      dos(sprintf('chmod +x %s', exefile));
      disp('SimBio is calculating the LeadFields, this may take some time ...')
      
      stopwatch = tic;
      dos(['./' exefile]);
      disp([ 'elapsed time: ' num2str(toc(stopwatch)) ])
      
      [lf] = sb_read_msr(outfile);
      cleaner(exefile,elcfile,dipfile,parfile,meshfile,outfile)
      cd(tmpfolder)
    end
  end
  
catch
  warning('an error occurred while running SimBio');
  rethrow(lasterror)
  cleaner(dipfile,elcfile,outfile,exefile)
  cd(tmpfolder)
end

function cleaner(exefile,elcfile,dipfile,parfile,meshfile,outfile)
delete(exefile);
delete(elcfile);
delete(dipfile);
delete(parfile);
delete(meshfile);
delete(outfile);
delete([outfile '.fld']);
delete([outfile '.hex']);
delete([outfile '.pot']);
delete([outfile '.pts']);
delete(['*.v.ascii']);
delete(['*.v.potential']);
