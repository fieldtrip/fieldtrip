function vol = ft_headmodel_fem_simbio(seg,varargin)
% FT_HEADMODEL_FEM_SIMBIO reads a volume conduction model from a Vista .v
% file
%
%
% Use as
%   vol = ft_headmodel_fem_simbio(seg)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD, TEST_HEADMODEL_SIMBIO

% Copyright (C) 2011, Cristiano Micheli
ft_hastoolbox('simbio', 1);

wfmethod     = ft_getopt(varargin, 'wfmethod', 'cubes'); % wireframe method (default: cubes)
transform    = ft_getopt(varargin, 'transform', []);  % contains the tr. matrix from voxels to head coordinates
unit         = ft_getopt(varargin, 'unit', 'mm');     % contains the units of the transformation matrix (mm, cm, ...)
tissue       = ft_getopt(varargin, 'tissue', []);     % contains the labels of the tissues
tissueval    = ft_getopt(varargin, 'tissueval', []);  % contains the tissue values (an integer for each compartment)
tissuecond   = ft_getopt(varargin, 'tissuecond', []); % contains the tissue conductivities
tissueweight = ft_getopt(varargin, 'tissueweight', ones(1,numel(tissueval))); % contains the weigths for vgrid (how dense the wireframe of each tissue conpartment should be done)
deepelec     = istrue(ft_getopt(varargin, 'deepelec', 'no'));% used in the case of deep voxel solution
% bnd          = ft_getopt(varargin, 'bnd', []);        % used in the case the solution has to be calculated on a triangulated surface (typically the scalp)
sens         = ft_getopt(varargin, 'sens', []);

% define a wireframe for each compartment (FEM grid)
wf = []; tMat = [];

if isempty(sens)
  error('A set of sensors is required')
end

if strcmp(wfmethod,'cubes')
  if isunix
    try % write wireframe grid with vgrid
      tmpfolder = cd;
      cd(tempdir)
      
      % temporary file names for vgrid call
      [tmp,tname] = fileparts(tempname);
      shfile    = [tname '.sh'];
      [tmp,tname] = fileparts(tempname);
      meshfile  = [tname '.v'];
      [tmp,tname] = fileparts(tempname);
      MRfile = [tname '.v'];
      [tmp,tname] = fileparts(tempname);
      materialsfile = [tname '.mtr'];
      % temporary file names for simbio call
      [tmp,tname] = fileparts(tempname);
      exefile   = [tname '.sh'];
      [tmp,tname] = fileparts(tempname);
      transfermatrix = [tname];
      [tmp,tname] = fileparts(tempname);
      parfile   = [tname '.par'];
      [tmp,tname] = fileparts(tempname);
      elcfile   = [tname '.elc'];
      
      if ~ft_hastoolbox('simbio',1,0) || ~ft_hastoolbox('vgrid',1,0)
        error('Cannot proceed without Simbio and Vgrid toolboxes')
      end
      
      % write the segmented volume in a Vista format .v file
      write_vista_vol(size(seg), seg, MRfile);
      
      % write the materials file (assign tissue values to the elements of the FEM grid)
      % see tutorial http://www.rheinahrcampus.de/~medsim/vgrid/manual.html
      sb_write_materials(materialsfile,tissueval,tissue,tissueweight);
      
      % write the shell file
      efid  = fopen(shfile, 'w');
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['vgrid -in ' MRfile ' -out ' meshfile ' -min 1 -max 1 -elem cube', ...
        ' -material ' materialsfile ' -smooth shift -shift 0.49 2>&1 > /dev/null\n']);
      fclose(efid);
      dos(sprintf('chmod +x %s', shfile));
      disp('Vgrid is writing the wireframe mesh file, this may take some time ...')
      stopwatch = tic;
      
      % Use vgrid to get the wireframe
      dos(['./' shfile]);
      disp([ 'elapsed time: ' num2str(toc(stopwatch)) ])
      % FIXME: think about adding a translation due to conversion between indices and vertices' world coordinates
      
      % calculate the FE transfer matrix
      
      % read the mesh points and assign wireframe (also called 'FEM grid', or '3d mesh')
      if ft_hastoolbox('fileio',1)
        [wf] = ft_read_headshape(meshfile);
      else
        error('you need the fileio module!\n')
      end
      
      if deepelec
        sb_write_elc(warp_apply(inv(transform),sens.chanpos),sens.label,elcfile,1);
      else
        sb_write_elc(warp_apply(inv(transform),sens.chanpos),sens.label,elcfile);
      end
      
      % write parfile
      disp('writing the parameters file on disk...')
      sb_write_par(parfile,'cond',tissuecond,'labels',tissueval);
      
      % write exefile
      efid = fopen(exefile, 'w');
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['ipm_linux_opt -i FEtransfermatrix -h ./' meshfile ' -s ./' elcfile, ...
        ' -o ./' transfermatrix ' -p ./' parfile ' -sens EEG \n']); %2>&1 > /dev/null
      fclose(efid);
      dos(sprintf('chmod +x %s', exefile));
      disp('simbio is calculating the transfer matrix, this may take some time ...')
      
      % generate the transfer matrix
      stopwatch = tic;
      dos(['./' exefile]);
      disp([ 'elapsed time: ' num2str(toc(stopwatch)) ])
      
      % read the transfer matrix
      transfer = sb_read_transfer(transfermatrix);
      cleaner(shfile,meshfile,MRfile,materialsfile,exefile,transfermatrix,parfile,elcfile)
    catch ME
      disp('The transfer matrix was not written')
      cleaner(shfile,meshfile,MRfile,materialsfile,exefile,transfermatrix,parfile,elcfile)
      cd(tmpfolder)
      rethrow(ME)
    end
    
  end %  if isunix
  
elseif strcmp(wfmethod,'tetra')
  % not yet implemented
else
  error('Unsupported method')
end

vol = [];
vol.wf        = wf;
vol.cond      = tissuecond;
vol.transform = transform;
vol.unit      = unit;
vol.type      = 'simbio';
vol.transfer  = transfer;

% if ~isempty(bnd)
%   vol.bnd        = bnd;
% else
if ~isempty(deepelec)
  vol.deepelec  = deepelec;
end

function cleaner(shfile,meshfile,MRfile,materialsfile,exefile,transfermatrix,parfile,elcfile)
delete(shfile);
delete(meshfile);
delete(MRfile);
delete(materialsfile);
delete(exefile);
delete([transfermatrix '.bin']);
delete([transfermatrix '.mat']);
delete(parfile);
delete(elcfile);
