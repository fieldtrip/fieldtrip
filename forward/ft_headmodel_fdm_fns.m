function vol = ft_headmodel_fdm_fns(varargin)

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
segmentation = ft_getopt(varargin, 'segmentation', []);
tissue       = ft_getopt(varargin, 'tissue', []);
tissueval    = ft_getopt(varargin, 'tissueval', []);
tissuecond   = ft_getopt(varargin, 'tissuecond', []);

% load the default cond matrix in case not specified
if isempty(condmatrix)
  if max(tissueval)<=12
    condmatrix = fns_contable_write;
  else
    % all tissue not listed as defaults (This is a bit rigid!!)
    sel = find(tissueval>12);
    newtissue     = tissue(sel);
    newtissueval  = tissueval(sel);
    newtissuecond = tissuecond(sel);    
    condmatrix = fns_contable_write('newtissue',newtissue,'newtissueval',newtissueval,'newtissuecond',newtissuecond);
  end
end

% check the consistency between tissue values and the segmentation
vecval = ismember(tissueval,unique(segmentation(:)));
if any(vecval)==0
  warning('Some of the tissue values are not in the segmentation')
end

% start with an empty volume conductor
vol = [];
vol.condmatrix = condmatrix;
vol.seg        = segmentation; 
vol.tissue     = tissue;
vol.tissueval  = tissueval;
vol.type       = 'fns';

