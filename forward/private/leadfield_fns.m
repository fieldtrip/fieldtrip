function [lf] = leadfield_fns(dip, elc, vol)

% LEADFIELD_FNS leadfields for a set of dipoles
%
% [lf] = leadfield_fns(dip, elc, vol);
%
% with input arguments
%   dip     positions of the dipoles (matrix of dimensions NX3)
%   elc     positions of the electrodes (matrix of dimensions MX3)
% 
% and vol being a structure with the element
%   vol.meshfile file containing the 3D mesh filename
%
% The output lf is the leadfields matrix of dimensions M (rows) X N*3 (cols)

% Copyright (C) 2011, Cristiano Micheli and Hung Dang

% store the current path and change folder to the temporary one
tmpfolder = cd;

if ft_senstype(elc, 'meg')
  error('FNS solver works for EEG data only!')
end

try
  if ~ispc
    cd(tempdir)
    [junk,tname] = fileparts(tempname);
    exefile = [tname '.sh'];   
    [junk,tname] = fileparts(tempname);
    elecfile   = [tname '.h5'];
    [junk,tname] = fileparts(tempname);
    confile   = [tname '.csv'];
    [junk,tname] = fileparts(tempname);
    datafile  = [tname '.h5'];
    
    tolerance = 1e-8; % FIXME: initialize as input argument

    % write the electrodes positions and the conductivity file in the temporary folder
    fns_elec_write(elc.pnt,vsize,dimpos,elecfile); % FIXME: this does not work yet
    fns_contable_write(vol.cond,confile);
    % fns_write_dip(dip,dipfile); % FNS ideally this should work like this!
    
    % Exe file
    % FIXME: does FNS have a switch for parallel processing (to run on more cores)?
    efid = fopen(exefile, 'w');
    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['elecsfwd --img ' vol.segfile ' --electrodes ./' elecfile ' --data ./', ...
                   datafile ' --contable ./' confile ' --TOL ' num2str(tolerance) ' 2>&1 > /dev/null\n']);
    end
    fclose(efid);
    
    dos(sprintf('chmod +x %s', exefile));
    dos(['./' exefile]);
    
    cleaner(elecfile,datafile,confile,exefile)
    
    % FNS calculates all solutions in all nodes (=voxels) at a time
    % the following operations extract the system matrix solution for gray
    % matter (FIXME: this code needs revision from FNS party)
    [junk,tname] = fileparts(tempname);
    regfile      = [tname '.h5'];
    exestr = sprintf('./img_get_gray --img %s --rgn %s',vol.segfile,regfile);
    system(exestr)
    
    % The following code extracts the matrix and uses it to calculate the
    % leadfields
    % FIXME: this code needs revision
    [junk,tname] = fileparts(tempname);
    recipfile    = [tname '.h5'];
    exestr = sprintf('./extract_data --data %s --rgn %s --newdata %s',datafile,regfile,recipfile);
    system(exestr)
    [data,compress,gridlocs,node_sizes,voxel_sizes] = fns_read_recipdata(recipfile);
    [lf] = fns_leadfield(data,compress,node_sizes,dip);

    cd(tmpfolder)
  end

catch
  warning('an error occurred while running FNS');
  rethrow(lasterror)
  cleaner(dipfile,elcfile,outfile,exefile)
  cd(tmpfolder)
end

function cleaner(elecfile,datafile,confile,exefile,regfile,recipfile)
delete(elecfile);
delete(datafile);
delete(confile);
delete(exefile);
delete(regfile);
delete(recipfile);
