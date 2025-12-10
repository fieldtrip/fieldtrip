function nii = niiheader2jnii(nii0)

nii = struct();
nii.NIFTIHeader.NIIHeaderSize =  nii0.hdr.sizeof_hdr;
if (isfield(nii0.hdr, 'data_type'))
    nii.NIFTIHeader.A75DataTypeName = deblank(char(nii0.hdr.data_type));
    nii.NIFTIHeader.A75DBName =      deblank(char(nii0.hdr.db_name));
    nii.NIFTIHeader.A75Extends =     nii0.hdr.extents;
    nii.NIFTIHeader.A75SessionError = nii0.hdr.session_error;
    nii.NIFTIHeader.A75Regular =     nii0.hdr.regular;
end
nii.NIFTIHeader.DimInfo.Freq =   bitand(uint8(nii0.hdr.dim_info), 7);
nii.NIFTIHeader.DimInfo.Phase =  bitand(bitshift(uint8(nii0.hdr.dim_info), -3), 7);
nii.NIFTIHeader.DimInfo.Slice =  bitand(bitshift(uint8(nii0.hdr.dim_info), -6), 7);
nii.NIFTIHeader.Dim =            nii0.hdr.dim(2:2 + nii0.hdr.dim(1) - 1);
nii.NIFTIHeader.Param1 =         nii0.hdr.intent_p1;
nii.NIFTIHeader.Param2 =         nii0.hdr.intent_p2;
nii.NIFTIHeader.Param3 =         nii0.hdr.intent_p3;
nii.NIFTIHeader.Intent =         niicodemap('intent', nii0.hdr.intent_code);
nii.NIFTIHeader.DataType =       niicodemap('datatype', nii0.hdr.datatype);
nii.NIFTIHeader.BitDepth =       nii0.hdr.bitpix;
nii.NIFTIHeader.FirstSliceID =   nii0.hdr.slice_start;
nii.NIFTIHeader.VoxelSize =      nii0.hdr.pixdim(2:2 + nii0.hdr.dim(1) - 1);
nii.NIFTIHeader.Orientation =    struct('x', 'r', 'y', 'a', 'z', 's');
if (nii0.hdr.pixdim(1) < 0)
    nii.NIFTIHeader.Orientation =    struct('x', 'l', 'y', 'a', 'z', 's');
end
nii.NIFTIHeader.NIIByteOffset =  nii0.hdr.vox_offset;
nii.NIFTIHeader.ScaleSlope =     nii0.hdr.scl_slope;
nii.NIFTIHeader.ScaleOffset =    nii0.hdr.scl_inter;
nii.NIFTIHeader.LastSliceID =    nii0.hdr.slice_end;
nii.NIFTIHeader.SliceType =      niicodemap('slicetype', nii0.hdr.slice_code);
nii.NIFTIHeader.Unit.L =         niicodemap('unit', bitand(uint8(nii0.hdr.xyzt_units), 7));
nii.NIFTIHeader.Unit.T =         niicodemap('unit', bitand(uint8(nii0.hdr.xyzt_units), 56));
nii.NIFTIHeader.MaxIntensity =   nii0.hdr.cal_max;
nii.NIFTIHeader.MinIntensity =   nii0.hdr.cal_min;
nii.NIFTIHeader.SliceTime =      nii0.hdr.slice_duration;
nii.NIFTIHeader.TimeOffset =     nii0.hdr.toffset;
if (isfield(nii0.hdr, 'glmax'))
    nii.NIFTIHeader.A75GlobalMax =   nii0.hdr.glmax;
    nii.NIFTIHeader.A75GlobalMin =   nii0.hdr.glmin;
end
nii.NIFTIHeader.Description =    deblank(char(nii0.hdr.descrip));
nii.NIFTIHeader.AuxFile =        deblank(char(nii0.hdr.aux_file));
nii.NIFTIHeader.QForm =          nii0.hdr.qform_code;
nii.NIFTIHeader.SForm =          nii0.hdr.sform_code;
nii.NIFTIHeader.Quatern.b =      nii0.hdr.quatern_b;
nii.NIFTIHeader.Quatern.c =      nii0.hdr.quatern_c;
nii.NIFTIHeader.Quatern.d =      nii0.hdr.quatern_d;
nii.NIFTIHeader.QuaternOffset.x = nii0.hdr.qoffset_x;
nii.NIFTIHeader.QuaternOffset.y = nii0.hdr.qoffset_y;
nii.NIFTIHeader.QuaternOffset.z = nii0.hdr.qoffset_z;
nii.NIFTIHeader.Affine(1, :) =    nii0.hdr.srow_x;
nii.NIFTIHeader.Affine(2, :) =    nii0.hdr.srow_y;
nii.NIFTIHeader.Affine(3, :) =    nii0.hdr.srow_z;
nii.NIFTIHeader.Name =           deblank(char(nii0.hdr.intent_name));
nii.NIFTIHeader.NIIFormat =      deblank(char(nii0.hdr.magic));
if (isfield(nii0.hdr, 'extension'))
    nii.NIFTIHeader.NIIExtender =    nii0.hdr.extension;
end
nii.NIFTIHeader.NIIQfac_ =       nii0.hdr.pixdim(1);
nii.NIFTIHeader.NIIEndian_ =     nii0.endian;
if (isfield(nii0.hdr, 'reserved'))
    nii.NIFTIHeader.NIIUnused_ = nii0.hdr.reserved;
end
