function headmodel = ft_headmodel_fns(seg, varargin)

% FT_HEADMODEL_FNS creates the volume conduction structure to be used
% in the FNS forward solver.
%
% Use as
%   headmodel = ft_headmodel_fns(seg, ...)
%
% Optional input arguments should be specified in key-value pairs and
% can include
%   tissuecond       = matrix C [9XN tissue types]; where N is the number of
%                      tissues and a 3x3 tensor conductivity matrix is stored
%                      in each column.
%   tissue           = see fns_contable_write
%   tissueval        = match tissues of segmentation input
%   transform        = 4x4 transformation matrix (default eye(4))
%   sens             = sensor information (for which ft_datatype(sens,'sens')==1)
%   deepelec         = used in the case of deep voxel solution
%   tolerance        = scalar (default 1e-8)
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
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

ft_hastoolbox('fns', 1);

% get the optional arguments
tissue       = ft_getopt(varargin, 'tissue', []);
tissueval    = ft_getopt(varargin, 'tissueval', []);
tissuecond   = ft_getopt(varargin, 'tissuecond', []);
transform    = ft_getopt(varargin, 'transform', eye(4));
unit         = ft_getopt(varargin, 'unit', 'mm');
sens         = ft_getopt(varargin, 'sens', []);
deepelec     = ft_getopt(varargin, 'deepelec', []); % used in the case of deep voxel solution
tolerance    = ft_getopt(varargin, 'tolerance', 1e-8);

if isempty(sens)
  ft_error('A set of sensors is required')
end

if ispc
  ft_error('FNS only works on Linux and OS X')
end

% check the consistency between tissue values and the segmentation
vecval = ismember(tissueval,unique(seg(:)));
if any(vecval)==0
  ft_warning('Some of the tissue values are not in the segmentation')
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
  
  % this requires the fieldtrip/fileio toolbox
  ft_hastoolbox('fileio', 1);
  
  % create a fake mri structure and write the segmentation on disk
  disp('writing the segmentation file...')
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
  pos = ft_warp_apply(inv(transform),sens.elecpos); % in voxel coordinates!
  
  % convert pos into int32 datatype.
  hdf5write(elecfile, '/electrodes/gridlocs', int32(pos));
  
  % Exe file
  efid = fopen(exefile, 'w');
  if ~ispc
    fprintf(efid,'#!/usr/bin/env bash\n');
    fprintf(efid,['elecsfwd1 -img ' segfile ' -electrodes ./' elecfile ' -data ./', ...
      datafile ' -contable ./' confile ' -TOL ' num2str(tolerance) ' \n']); %2>&1 > /dev/null
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
headmodel = [];
headmodel.tissue     = tissue;
headmodel.tissueval  = tissueval;
headmodel.transform  = transform;
headmodel.unit       = unit;
headmodel.segdim     = size(seg);
headmodel.type       = 'fns';
headmodel.transfer   = transfer;

if ~isempty(deepelec)
  headmodel.deepelec  = deepelec;
end

function cleaner(segfile,confile,elecfile,exefile,datafile)
delete([segfile '.hdr']);
delete([segfile '.img']);
delete(confile);
delete(elecfile);
delete(exefile);
delete(datafile);
