function varargout = loadnifti (varargin)
%
%    jnii=loadnifti(filename)
%        or
%    nii=loadnifti(filename,option)
%
%    Read a NIfTI-1/2 (*.nii/.nii.gz) or Analyze 7.5 (*.hdr/*.img/.hdr.gz/.img.gz) 
%    image file.
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    Please run `help nii2jnii` to see the input output outputs.
%    This function is an alias to nii2jnii
%
%
%    this file is part of JNIfTI specification: https://github.com/fangq/jnifti
%
%    License: Apache 2.0, see https://github.com/fangq/jnifti for details
%

[varargout{1:nargout}]=nii2jnii(varargin{:});