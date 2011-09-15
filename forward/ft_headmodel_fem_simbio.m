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
unit         = ft_getopt(varargin, 'unit', 'cm');     % contains the units of the transformation matrix (mm, cm, ...)
tissue       = ft_getopt(varargin, 'tissue', []);     % contains the labels of the tissues
tissueval    = ft_getopt(varargin, 'tissueval', []);  % contains the tissue values (an integer for each compartment)
tissuecond   = ft_getopt(varargin, 'tissuecond', []); % contains the tissue conductivities
tissueweight = ft_getopt(varargin, 'tissueweight', ones(1,numel(tissueval))); % contains the weigths for vgrid (how dense the wireframe of each tissue conpartment should be done)
deepelec     = ft_getopt(varargin, 'deepelec', []);   % used in the case of deep voxel solution
bnd          = ft_getopt(varargin, 'bnd', []);        % used in the case the solution has to be calculated on a triangulated surface (typically the scalp)

% define a wireframe for each compartment (FEM grid)
wf = []; tMat = [];
if strcmp(wfmethod,'cubes')
  tmpfolder = cd;
  cd(tempdir)
  % write the segmented volume in a Vista format .v file
  [~,tname] = fileparts(tempname);
  MRfile = [tname '.v'];
  % write a vista volumetric file
  %   ft_write_volume(MRfile, seg,'dataformat', 'vista');
  write_vista_vol(size(seg), seg, MRfile);
  
  if ~ft_hastoolbox('simbio')
    error('Cannot write a materials file without the Vista/Simbio toolbox')
  end
  % write the materials file (assign tissue values to the elements of the FEM grid)
  % see tutorial: http://www.rheinahrcampus.de/~medsim/vgrid/manual.html
  [~,tname] = fileparts(tempname);
  materialsfile = [tname '.mtr'];
  sb_write_materials(materialsfile,tissueval,tissue,tissueweight);

  % Use vgrid to get the wireframe
  %   vroot = '/home/coherence/crimic/test/SimBio/fromJohannes/vgrid1.3.1/program/';
  [~,tname] = fileparts(tempname);
  meshfile  = [tname '.v'];
  [~,tname] = fileparts(tempname);
  shfile    = [tname '.sh'];
  efid  = fopen(shfile, 'w');
  fprintf(efid,'#!/usr/bin/env bash\n');
  % works with voxels
  fprintf(efid,['vgrid -in ' MRfile ' -out ' meshfile ' -min 1 -max 1 -elem cube', ...
    ' -material ' materialsfile ' -smooth shift -shift 0.49 2>&1 > /dev/null\n']);
  fclose(efid);
  % execute command
  dos(sprintf('chmod +x %s', shfile));
  disp('Vgrid is writing the wireframe mesh file, this may take some time ...')
  stopwatch = tic;
  % FIXME: think about: add a translation due to conversion between indices and vertices'
  % worlds coordinates?
  
  try
    dos(['./' shfile]);
    disp([ 'elapsed time: ' num2str(toc(stopwatch)) ])
  catch
    disp('Error in writing the wireframe mesh file')
    cleaner(MRfile,materialsfile,meshfile,shfile)
    rethrow(ME)
  end
  
  try
    % read the mesh points
    [nodes,elements,labels] = read_vista_mesh(meshfile);
    % assign wireframe (also called 'FEM grid', or '3d mesh')
    wf.nd = nodes;
    wf.el = elements;
    wf.labels = labels;
    
    % exe file
    if ~ispc
      [~,tname] = fileparts(tempname);
      exefile   = [tname '.sh'];
      [~,tname] = fileparts(tempname);
      transfermatrix = [tname];
      [~,tname] = fileparts(tempname);
      parfile   = [tname '.par'];
      
      % write parfile
      disp('writing the parameters file on disk...')
      sb_write_par(parfile,'cond',tissuecond,'labels',tissueval); %unique(vol.wf.labels)
      
      efid = fopen(exefile, 'w');
      fprintf(efid,'#!/usr/bin/env bash\n');
      % generate the transfer matrix 
      fprintf(efid,['ipm_linux_opt -i FEtransfermatrix -h ./' meshfile ' -s ./' meshfile, ...
        ' -o ./' transfermatrix ' -p ./' parfile ' -sens EEG 2>&1 > /dev/null\n']);
      fclose(efid);
      
      dos(sprintf('chmod +x %s', exefile));
      disp('simbio is calculating the transfer matrix, this may take some time ...')

      stopwatch = tic;
      dos(['./' exefile]);
      disp([ 'elapsed time: ' num2str(toc(stopwatch)) ])
      % read the tr matrix 
      tMat = sb_read_transfer(transfermatrix);
      tMat = tMat.mat;
    end
    
    cleaner(MRfile,materialsfile,meshfile,shfile,exefile,parfile,transfermatrix) 
  catch
    disp('mesh was not written, check the presence of vgrid in the path')
    cleaner(MRfile,materialsfile,meshfile,shfile,exefile,parfile,transfermatrix)
    rethrow(ME)
  end
  cd(tmpfolder)
  
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
vol.mat       = tMat;

if ~isempty(bnd)
  vol.bnd        = bnd;
elseif ~isempty(deepelec)
  vol.deepelec  = deepelec;
end

function cleaner(MRfile,materialsfile,meshfile,shfile,exefile,parfile,transfermatrix)
delete(MRfile);
delete(materialsfile);
delete(meshfile);
delete(shfile);
delete(exefile);
delete(parfile);
delete(transfermatrix);
