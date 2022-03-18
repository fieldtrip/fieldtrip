function savebnii(jnii, filename, varargin)
%
%    savebnii(jniidata, bniifile)
%       or
%    savebnii(jniidata, bniifile, 'Param1',value1, 'Param2',value2,...)
%
%    Save an in-memory JNIfTI structure into a binary-JNIfTI file with format
%    defined in JNIfTI specification: https://github.com/fangq/jnifti
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        jnii: a structure (array) or cell (array). The data structure can
%            be completely generic or auxilary data without any JNIfTI
%            constructs. However, if a JNIfTI object is included, it shall
%            contain the below subfields (can appear within any depth of the
%            structure)
%                jnii.NIFTIHeader -  a structure containing the 1-to-1 mapped NIFTI-1/2 header
%                jnii.NIFTIData - the main image data array
%                jnii.NIFTIExtension - a cell array contaiing the extension data buffers
%        filename: the output file name to the binary-JNIfTI file (.bnii)
%        options: (optional) if saving to .bnii, please see the 
%               supported options for savebj.m (part of JSONLab).
%
%    example:
%        jnii=jnifticreate(uint8(magic(10)),'Name','10x10 magic matrix');
%        savebnii(jnii, 'magic10.bnii')
%        savebnii(jnii, 'magic10_debug.bnii','Debug',1)
%
%    this file is part of JNIfTI specification: https://github.com/fangq/jnifti
%
%    License: Apache 2.0, see https://github.com/fangq/jnifti for details
%

if(nargin<2)
    error('you must provide data and output file name');
end

if(~exist('savebj','file'))
    error('you must first install JSONLab from http://github.com/fangq/jsonlab/');
end

savebj('',jnii,'FileName',filename,varargin{:});