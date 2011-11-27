function vol = ft_headmodel_fdm_fns(seg,varargin)

% FT_HEADMODEL_FDM_FNS creates the volume conduction structure to be used 
% in the FNS forward solver.
%
% Use as
%   vol = ft_headmodel_fdm_fns(varargin)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = matrix C [9XN tissue types]
%
% The conductivity matrix C is stored in a 9xN table, where N is the number of tissues and a 
% 3x3 tensor conductivity matrix is stored in each row in column oriented format
% 
% Standard default values for conductivity matrix C are derived from 
% Saleheen HI, Ng KT. New finite difference formulations for general
% inhomogeneous anisotropic bioelectric problems. IEEE Trans Biomed Eng.
% 1997
% 
% Additional documentation available at:
% http://hunghienvn.nmsu.edu/wiki/index.php/FNS
% 
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2011, Cristiano Micheli and Hung Dang

ft_hastoolbox('fns', 1);

% get the optional arguments
tissue       = ft_getopt(varargin, 'tissue', []);
tissueval    = ft_getopt(varargin, 'tissueval', []);
tissuecond   = ft_getopt(varargin, 'tissuecond', []);
transform    = ft_getopt(varargin, 'transform', eye(4));
units        = ft_getopt(varargin, 'units', 'cm');
sens         = ft_getopt(varargin, 'sens', []);
deepelec     = ft_getopt(varargin, 'deepelec', []); % used in the case of deep voxel solution
tolerance    = ft_getopt(varargin, 'tolerance', 1e-8);

if isempty(sens)
  error('A set of sensors is required')
end

if ispc
  error('FNS only works on Linux and OSX')
end

% check the consistency between tissue values and the segmentation
vecval = ismember(tissueval,unique(seg(:)));
if any(vecval)==0
  warning('Some of the tissue values are not in the segmentation')
end

% create the files to be written
try
    tmpfolder = pwd;
    
    cd(tempdir)
    [tmp,tname] = fileparts(tempname);
    segfile   = [tname];     
    [tmp,tname] = fileparts(tempname);
    confile   = [tname '.csv'];
    [tmp,tname] = fileparts(tempname);
    elecfile = [tname '.h5']; 
    [tmp,tname] = fileparts(tempname);
    exefile   = [tname '.sh'];     
    [tmp,tname] = fileparts(tempname);
    datafile  = [tname '.h5'];
       
    % create a fake mri structure and write the segmentation on disk
    disp('writing the segmentation file...')
    if ~ft_hastoolbox('fileio')
      error('You must have the fileio module to go on')
    end    
    mri = [];
    mri.dim = size(seg);
    mri.transform = eye(4);
    mri.seg = uint8(seg);
    
    cfg = [];
    cfg.datatype = 'uint8';
    cfg.coordsys  = 'ctf';
    cfg.parameter = 'seg';
    cfg.filename  = segfile;
    cfg.filetype  = 'analyze';
    ft_volumewrite(cfg, mri);     
    
    % write the cond matrix on disk, load the default cond matrix in case not specified
    disp('writing the conductivity file...')
    condmatrix = fns_contable_write('tissue',tissue,'tissueval',tissueval,'tissuecond',tissuecond);
    csvwrite(confile,condmatrix);

    % write the positions of the electrodes on disk
    disp('writing the electrodes file...')
    pos = warp_apply(inv(transform),sens.chanpos); % in voxel coordinates!

    % convert pos into int32 datatype. 
    hdf5write(elecfile, '/electrodes/gridlocs', int32(pos));
    
    % Exe file 
    efid = fopen(exefile, 'w');
    if ~ispc
      fprintf(efid,'#!/usr/bin/env bash\n');
      fprintf(efid,['elecsfwd1 -img ' segfile ' -electrodes ./' elecfile ' -data ./', ...
                   datafile ' -contable ./' confile ' -TOL ' num2str(tolerance) ' \n']);%2>&1 > /dev/null
    end
    fclose(efid);
    
    % run the shell instructions
    dos(sprintf('chmod +x %s', exefile));
    dos(['./' exefile]);
    
    % FIXME: find a cleverer way to store the huge transfer matrix (vista?)
    [transfer,status] = fns_read_transfer(datafile);
    
    cleaner(segfile,confile,elecfile,exefile,datafile)
    
catch ME
    disp('The transfer matrix was not written')
    cleaner(segfile,confile,elecfile,exefile,datafile)
    cd(tmpfolder)
    rethrow(ME)
end

% start with an empty volume conductor
vol = [];
vol.tissue     = tissue;
vol.tissueval  = tissueval;
vol.transform  = transform;
vol.segdim     = size(seg);
vol.units      = units;
vol.type       = 'fns';
vol.transfer   = transfer;

if ~isempty(deepelec)
  vol.deepelec  = deepelec;
end

function cleaner(segfile,confile,elecfile,exefile,datafile)
delete([segfile '.hdr']);
delete([segfile '.img']);
delete(confile);
delete(elecfile);
delete(exefile);
delete(datafile);
