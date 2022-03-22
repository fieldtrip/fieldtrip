function savejnifti(jnii, filename, varargin)
%
%    savejnifti(jnii, outputfile)
%       or
%    savejnifti(jnii, outputfile, 'Param1',value1, 'Param2',value2,...)
%
%    Save an in-memory JNIfTI structure into a JNIfTI file with format
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
%        outputfile: the output file name to the JNIfTI file 
%                *.bnii for binary JNIfTI file
%                *.jnii for text JNIfTI file
%        options: (optional) if saving to a .bnii file, please see the options for
%               savebj.m (part of JSONLab); if saving to .jnii, please see the 
%               supported options for savejson.m (part of JSONLab).
%
%    example:
%        jnii=jnifticreate(uint8(magic(10)),'Name','10x10 magic matrix');
%        savejnifti(jnii, 'magic10.jnii')
%        savejnifti(jnii, 'magic10_debug.bnii','Debug',1)
%
%    this file is part of JNIfTI specification: https://github.com/fangq/jnifti
%
%    License: Apache 2.0, see https://github.com/fangq/jnifti for details
%

if(nargin<2)
    error('you must provide data and output file name');
end

if(~exist('savejson','file'))
    error('you must first install JSONLab from http://github.com/fangq/jsonlab/');
end

if(regexp(filename,'\.jnii$'))
    savejnii(jnii,filename,varargin{:});
elseif(regexp(filename,'\.bnii$'))
    savebnii(jnii,filename,varargin{:});
else
    error('file suffix must be .jnii for text JNIfTI or .bnii for binary JNIfTI');
end
