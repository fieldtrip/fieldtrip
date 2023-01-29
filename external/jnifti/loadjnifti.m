function jnii=loadjnifti(filename, varargin)
%
%    jnii=loadjnifti(inputfile)
%       or
%    jnii=loadjnifti(inputfile, 'Param1',value1, 'Param2',value2,...)
%
%    Load a standard NIFTI-1/2 file or text or binary JNIfTI file with
%    format defined in JNIfTI specification: https://github.com/fangq/jnifti
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        inputfile: the output file name to the JNIfTI or NIFTI-1/2 file
%                *.bnii for binary JNIfTI file
%                *.jnii for text JNIfTI file
%                *.nii  for NIFTI-1/2 files
%        options: (optional) if loading from a .bnii file, please see the options for
%               loadbj.m (part of JSONLab); if loading from a .jnii, please see the 
%               supported options for loadjson.m (part of JSONLab).
%
%    output:
%        jnii: a structure (array) or cell (array). The data structure can
%            be completely generic or auxilary data without any JNIfTI
%            constructs. However, if a JNIfTI object is included, it shall
%            contain the below subfields (can appear within any depth of the
%            structure)
%                jnii.NIFTIHeader -  a structure containing the 1-to-1 mapped NIFTI-1/2 header
%                jnii.NIFTIData - the main image data array
%                jnii.NIFTIExtension - a cell array contaiing the extension data buffers
%
%    example:
%        jnii=jnifticreate(uint8(magic(10)),'Name','10x10 magic matrix');
%        savejnifti(jnii, 'magic10.jnii')
%        newjnii=loadjnifti('magic10.jnii');
%
%    this file is part of JNIfTI specification: https://github.com/fangq/jnifti
%
%    License: Apache 2.0, see https://github.com/fangq/jnifti for details
%

if(nargin<1)
    error('you must provide data and output file name');
end

if(~exist('savejson','file'))
    error('you must first install JSONLab from http://github.com/fangq/jsonlab/');
end

if(regexp(filename,'\.nii$'))
    jnii=nii2jnii(filename,'jnii');
elseif(regexp(filename,'\.jnii$'))
    jnii=loadjson(filename,varargin{:});
elseif(regexp(filename,'\.bnii$'))
    jnii=loadbj(filename,varargin{:});
else
    error('file suffix must be .jnii for text JNIfTI, .bnii for binary JNIfTI or .nii for NIFTI-1/2 files');
end
