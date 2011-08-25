function [lf] = leadfield_fns(vx, vol, tol)

% LEADFIELD_FNS calculates the FDM forward solution for a set of voxels positions
%
% [lf] = leadfield_fns(vx, elc, vol);
%
% with input arguments
%   vx     positions of the output voxels (matrix of dimensions NX3)
%   vol    strucure of the volume conductor
%   tol    tolerance
% 
% The important elements of the vol structure are:
%   vol.condmatrix: a 9XT (T tissues) matrix containing the conductivities
%   vol.seg: a segmented/labelled MRI
%
% The output lf is the leadfields matrix of dimensions M (rows) X N voxels

% Copyright (C) 2011, Cristiano Micheli and Hung Dang

% store the current path and change folder to the temporary one
tmpfolder = cd;

try
  if ~ispc
    cd(tempdir)
    [junk,tname] = fileparts(tempname);
    exefile = [tname '.sh'];   
    [junk,tname] = fileparts(tempname);
    confile   = [tname '.csv'];
    [junk,tname] = fileparts(tempname);
    datafile  = [tname '.h5'];
    [junk,tname] = fileparts(tempname);
    segfile  = [tname];    
    [junk,tname] = fileparts(tempname);
    vxfile  = [tname '.csv'];
    
    if isempty(tol)
      tolerance = 1e-8;
    else
      tolerance = tol;
    end
    
    % creat a fake mri structure
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
    
    % write the positions of the voxels on disk
    csvwrite(vxfile,vx); 
    
    % Exe file (waiting for the parallel processing version of the binary)
    efid = fopen(exefile, 'w');
    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['elecsfwd --img ' segfile ' --electrodes ./' vxfile ' --data ./', ...
                   datafile ' --contable ./' confile ' --TOL ' num2str(tolerance) ' 2>&1 > /dev/null\n']);
    end
    fclose(efid);
    
    % run the shell instructions
    dos(sprintf('chmod +x %s', exefile));
    dos(['./' exefile]);
    cleaner(segfile,vxfile,datafile,confile,exefile)
    
% NOTE: this wont be necessary in the new version
%       % extract the system matrix solution for gray matter 
%       [junk,tname] = fileparts(tempname);
%       regfile      = [tname '.h5'];
%       [junk,tname] = fileparts(tempname);
%       recipfile    = [tname '.h5'];    
%       exestr = sprintf('./img_get_gray --img %s --rgn %s',vol.segfile,regfile);
%       system(exestr)
% 
%       % extract the matrix and uses it to calculate the leadfields
%       exestr = sprintf('./extract_data --data %s --rgn %s --newdata %s',datafile,regfile,recipfile);
%       system(exestr)
%       [data,compress,gridlocs,node_sizes,voxel_sizes] = fns_read_recipdata(recipfile);
%       [lf] = fns_leadfield(data,compress,node_sizes,dip);

    cd(tmpfolder)
  end

catch
  warning('an error occurred while running FNS');
  cleaner(segfile,vxfile,datafile,confile,exefile)
  rethrow(lasterror)
  cd(tmpfolder)
end

function cleaner(segfile,vxfile,datafile,confile,exefile)
delete(segfile);
delete(vxfile);
delete(datafile);
delete(confile);
delete(exefile);
