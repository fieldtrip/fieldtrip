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

% get the optional arguments
condmatrix   = ft_getopt(varargin, 'condmatrix', []);
tissue       = ft_getopt(varargin, 'tissue', []);
tissueval    = ft_getopt(varargin, 'tissueval', []);
tissuecond   = ft_getopt(varargin, 'tissuecond', []);
bnd          = ft_getopt(varargin, 'bnd', []);
transform    = ft_getopt(varargin, 'transform', eye(4));
units        = ft_getopt(varargin, 'units', 'cm');
deepelec     = ft_getopt(varargin, 'deepelec', []); % used in the case of deep voxel solution

if isempty(deepvoxel) && isempty(bnd)
  error('Either a deep electrode or a boundary have to be given')
end

% load the default cond matrix in case not specified
condmatrix = fns_contable_write('tissue',tissue,'tissueval',tissueval,'tissuecond',tissuecond);

% check the consistency between tissue values and the segmentation
vecval = ismember(tissueval,unique(seg(:)));
if any(vecval)==0
  warning('Some of the tissue values are not in the segmentation')
end

% start with an empty volume conductor
vol = [];
vol.condmatrix = condmatrix;
vol.seg        = seg; 
vol.tissue     = tissue;
vol.tissueval  = tissueval;
vol.transform  = transform;
vol.units      = units;
vol.type       = 'fns';

% bnd has always the precedence
if ~isempty(bnd)
  vol.bnd        = bnd;
else
  vol.deepelec  = deepelec;
end
