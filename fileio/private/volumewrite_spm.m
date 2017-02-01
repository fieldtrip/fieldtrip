function [Va] = volumewrite_spm(filename, data, transform, spmversion)

% VOLUMEWRITE_SPM writes anatomical or functional MRI volume data to analyze or nifti format
% using the SPM toolbox.
%
% Use as
%   [Va] = volumewrite_spm(filename, data, transform)

% Copyright (C) 2006, Robert Oostenveld
% Copyright (C) 2011, Jan-Mathijs Schoffelen

if nargin<4 || isempty(spmversion)
  spmversion = 'spm8';
end

% check whether the required SPM toolbox is available
ft_hastoolbox(spmversion, 1);

if 0
  % not all datatypes are supported by SPM, this can be checked with
  spm_type('double')  %           -- Double precision floating point number array
  spm_type('single')  %           -- Single precision floating point number array
  spm_type('logical') %           -- Logical array
  spm_type('char')    %           -- Character array
  spm_type('cell')    %           -- Cell array
  spm_type('struct')  %           -- Structure array
  spm_type('int8')    %           -- 8-bit signed integer array
  spm_type('uint8')   %           -- 8-bit unsigned integer array
  spm_type('int16')   %           -- 16-bit signed integer array
  spm_type('uint16')  %           -- 16-bit unsigned integer array
  spm_type('int32')   %           -- 32-bit signed integer array
  spm_type('uint32')  %           -- 32-bit unsigned integer array
  spm_type('int64')   %           -- 64-bit signed integer array
  spm_type('uint64')  %           -- 64-bit unsigned integer array
end

typ       = spm_type(class(data));
dim       = size(data);
if isnan(typ)
  % convert every unsupported data type into double
  data      = double(data);
  typ       = spm_type(class(data));
end

switch lower(spmversion)
  case 'spm2'
    %see spm_vol
    Va         = [];
    Va.mat     = transform;
    Va.fname   = filename;
    Va.dim     = [dim typ];
    Va.n       = 1;
    Va.pinfo   = [1 0 0]';
    Va.private.hdr.dime.datatype = typ;
    Va         = spm_create_vol(Va);
    Va         = spm_write_vol(Va,data);
    
  case 'spm8'
    Va         = [];
    Va.mat     = transform;
    Va.fname   = filename;
    if numel(dim)>3
      Va.dim = dim(1:3);
      Va.n   = dim(4:end);
    else
      Va.dim     = dim;
      Va.n       = 1;
    end
    Va.pinfo   = [1 0 0]';
    %Va.dt      = [typ 1]; % this is not necessary because assigned in spm_create_vol
    Va         = spm_create_vol(Va);
    Va         = spm_write_vol(Va,data);
    
  case 'spm12'
    N     = nifti;
    N.mat = transform;
    N.mat_intent = 'Aligned';
    N.dat = file_array(filename, dim, 'FLOAT32-LE');
    create(N);
    switch length(N.dat.dim)
      case 2
        N.dat(:,:)      = data;
      case 3
        N.dat(:,:,:)    = data;
      case 4
        N.dat(:,:,:, :) = data;
      otherwise
        error('Invalid output dimensions');
    end
    Va = spm_vol(N.dat.fname);
    
  otherwise
    error('unsupported SPM version requested');
end
