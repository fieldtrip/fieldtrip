function [lf] = leadfield_simbio(dip, elc, vol)

% leadfield_simbio leadfields for a set of dipoles
%
% [lf] = leadfield_simbio(elc, vol, headmeshfile);
%
% with input arguments
%   elc     positions of the electrodes (matrix of dimensions mx3)
%           contains the positions of the vertices on which
%           the potentials are calculated
%           there can be 'deep electrodes' too!
%   vol.wf  contains a wireframe structure (or fem grid) of 'elements'
%           which is needed to perform a fem analysis.
%   dip     a matrix of dipoles (of dimensions nx3)
%
% the output lf is the leadfields matrix of dimensions m (rows) x n*3 (cols)

% copyright (c) 2011, Cristiano Micheli

% store the current path and change folder to the temporary one
tmpfolder = cd;

% add the simbio folder to the path
if ~isunix
  error('this only works on linux and osx')
end

try 
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
  disp('writing the electrodes file on disk...\n')
  if ~isfield(vol,'deepelec')
    sb_write_elc(warp_apply(inv(vol.transform),elc.chanpos),elc.label,elcfile);
  else
    sb_write_elc(warp_apply(inv(vol.transform),elc.chanpos),elc.label,elcfile,1);
  end
  disp('writing the dipoles file on disk...\n')
  sb_write_dip(warp_apply(inv(vol.transform),dip),dipfile);
  
  % write the parameters file, contains tissues conductivities, the fe
  % grid and simbio call details, mixed together
  disp('writing the parameters file on disk...')
  sb_write_par(parfile,'cond',vol.cond,'labels',unique(vol.wf.index)-100);
  
  % write the vol.wf in a vista format .v file
  disp('writing the mesh file on disk...\n')
  if ft_hastoolbox('fileio',1)
    ft_write_headshape(meshfile,vol.wf,'format','vista');
  else
    error('you need the fileio module!\n')
  end
  
  % write the transfer matrix on disk
  disp('writing the transfer matrix on disk...')
  sb_write_transfer(vol.transfer,transfermatrix);
  
  if ~ispc
    efid = fopen(exefile, 'w');
    fprintf(efid,'#!/usr/bin/env bash\n');
    fprintf(efid,['ipm_linux_opt -i sourcesimulation  -h ./' meshfile ' -s ./' elcfile ' -t ./' transfermatrix, ...
      ' -dip ./' dipfile ' -o ./' outfile ' -p ./' parfile ' -fwd FEM -sens EEG 2>&1 > /dev/null\n']);
    fclose(efid);
    
%     fprintf(efid,['ipm_linux_opt_Venant -i sourcesimulation -h ./' meshfile ' -s ./' elcfile, ...
%       ' -dip ./' dipfile ' -o ./' outfile ' -p ./' parfile ' -fwd FEM -sens EEG 2>&1 > /dev/null\n']);

    dos(sprintf('chmod +x %s', exefile));
    disp('simbio is calculating the leadfields, this may take some time ...')

    stopwatch = tic;
    dos(['./' exefile]);
    disp([ 'elapsed time: ' num2str(toc(stopwatch)) ])

    % extract the lead-fields
    [lf] = sb_read_msr(outfile);
    cleaner(exefile,elcfile,dipfile,parfile,meshfile,outfile,transfermatrix)
    cd(tmpfolder)
  end
  
catch
  warning('an error occurred while running simbio');
  rethrow(lasterror)
  cleaner(exefile,elcfile,dipfile,parfile,meshfile,outfile,transfermatrix)
  cd(tmpfolder)
end

function cleaner(exefile,elcfile,dipfile,parfile,meshfile,outfile,transfermatrix)
delete(exefile);
delete(elcfile);
delete(dipfile);
delete(parfile);
delete(meshfile);
delete(outfile);
delete(transfermatrix);
delete([outfile '.fld']);
delete([outfile '.hex']);
delete([outfile '.pot']);
delete([outfile '.pts']);
delete(['*.v.ascii']);
delete(['*.v.potential']);
