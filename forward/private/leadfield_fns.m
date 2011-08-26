function [lf] = leadfield_fns(posin, vol, tol)

% LEADFIELD_FNS calculates the FDM forward solution for a set of given
% dipolar sources
%
% [lf] = leadfield_fns(posin, vol, tol);
%
% with input arguments
%   posin  positions of the voxels containing the dipolar sources (MX3 matrix)
%   vol    structure of the volume conductor
%   tol    tolerance
% 
% The output argument lf 
% 
% The key elements of the vol structure are:
%   vol.condmatrix a 9XT (T tissues) matrix containing the conductivities
%   vol.seg        a segmented/labelled MRI
%   vol.posout     positions of the output voxels (NX3 matrix)
%
% The output lf is the forward solution matrix and contains the leadfields, calculated in 
% the output voxels (vol.posout) positions (a NXM matrix)

% Copyright (C) 2011, Cristiano Micheli and Hung Dang

% store the current path and change folder to the temporary one
tmpfolder = cd;

try
  if ~ispc
    cd(tempdir)
    [~,tname] = fileparts(tempname);
    exefile   = [tname '.sh'];   
    [~,tname] = fileparts(tempname);
    confile   = [tname '.csv'];
    [~,tname] = fileparts(tempname);
    datafile  = [tname '.h5'];
    [~,tname] = fileparts(tempname);
    segfile   = [tname];    
    [~,tname] = fileparts(tempname);
    vxinfile  = [tname '.h5'];
    [~,tname] = fileparts(tempname);
    vxoutfile = [tname '.h5'];
    [~,tname] = fileparts(tempname);
    regfile   = [tname '.h5'];
    [~,tname] = fileparts(tempname);
    recipfile = [tname '.h5'];  
    
    if isempty(tol)
      tolerance = 1e-8;
    else
      tolerance = tol;
    end
    
    % create a fake mri structure
    mri = [];
    mri.dim = size(vol.seg);
    mri.transform = eye(4);
    mri.seg = vol.seg;
    % write the segmentation on disk
    cfg = [];
    cfg.coordsys  = 'ctf';
    cfg.parameter = 'seg';
    cfg.filename  = segfile;
    cfg.filetype  = 'analyze';
    ft_volumewrite(cfg, mri);

    % write the cond matrix on disk
    csvwrite(confile,vol.condmatrix);
    
    % write the positions of the output voxels on disk
    fns_elec_write(vol.posout,[1 1 1],size(mri.seg),vxoutfile); 
    
    % Exe file (the parallel processing version of the binary will soon be available)
    efid = fopen(exefile, 'w');
    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['elecsfwd --img ' segfile ' --electrodes ./' vxoutfile ' --data ./', ...
                   datafile ' --contable ./' confile ' --TOL ' num2str(tolerance) ' 2>&1 > /dev/null\n']);
    end
    fclose(efid);
    
    % run the shell instructions
    % NOTE: This requires FNS to be compiled/installed and the binaries folder to be already present in the 
    % linux PATH
    dos(sprintf('chmod +x %s', exefile));
    sprintf('Calculating the reciprocity matrix for %d nodes...\n (this may take some time)',size(vol.posout,1))
    dos(['./' exefile]);
    lf   = fns_leadfield(datafile,posin);
    
%     % extract the forward solution for the given dipoles (posin)
%     % write posin positions in XXX
%     exestr = sprintf('./img_getvoxels --img %s --rgn %s --arg XXX',segfile,regfile); % FIXME: fix this!
%     system(exestr)
%     exestr = sprintf('./extract_data --data %s --rgn %s --newdata %s',datafile,regfile,recipfile);
%     system(exestr)
%     lf   = fns_leadfield(recipfile,posin);

    cleaner(segfile,vxinfile,vxoutfile,datafile,confile,exefile,regfile,recipfile) 

    cd(tmpfolder)
  end

catch
  warning('an error occurred while running FNS');
  cleaner(segfile,vxinfile,vxoutfile,datafile,confile,exefile,regfile,recipfile) 
  rethrow(lasterror)
  cd(tmpfolder)
end

function cleaner(segfile,vxinfile,vxoutfile,datafile,confile,exefile,regfile,recipfile) 
delete(segfile);
delete(vxinfile);
delete(vxoutfile);
delete(datafile);
delete(confile);
delete(exefile);
delete(regfile);
delete(recipfile);
