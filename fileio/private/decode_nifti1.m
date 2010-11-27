function H = decode_nifti1(blob)
% function H = decode_nifti1(blob)
%
% Decodes a NIFTI-1 header given as raw 348 bytes (uint8) into a Matlab structure
% that matches the C struct defined in nifti1.h, with the only difference that the
% variable length arrays "dim" and "pixdim" are cut off to the right size, e.g., the
% "dim" entry will only contain the relevant elements: 
% dim[0..7]={3,64,64,18,x,x,x,x} in C would become dim=[64,64,18] in Matlab.
%
% WARNING: This function currently ignores endianness !!!

% (C) 2010 S.Klanke

if class(blob)~='uint8'
	error 'Bad type for blob'
end
if length(blob)~=348
	error 'Blob must be exactly 348 bytes long'
end

% see nift1.h for information on structure
H = [];

magic = char(blob(345:347));
if blob(348)~=0 | magic~='ni1' & magic~='n+1'
	error 'Not a NIFTI-1 header!';
end

H.sizeof_hdr = typecast(blob(1:4),'int32');
H.data_type  = cstr2matlab(blob(5:14));
H.db_name    = cstr2matlab(blob(15:32));
H.extents    = typecast(blob(33:36),'int32');
H.session_error = typecast(blob(37:38),'int16');
H.regular    = blob(39);
H.dim_info   = blob(40);

dim = typecast(blob(41:56),'int16');
H.dim = dim(2:dim(1)+1);
H.intent_p1   = typecast(blob(57:60),'single');
H.intent_p2   = typecast(blob(61:64),'single');
H.intent_p3   = typecast(blob(65:68),'single');
H.intent_code = typecast(blob(69:70),'int16');
H.datatype    = typecast(blob(71:72),'int16');
H.bitpix      = typecast(blob(73:74),'int16');
H.slice_start = typecast(blob(75:76),'int16');
pixdim = typecast(blob(77:108),'single');
H.qfac        = pixdim(1);
H.pixdim      = pixdim(2:dim(1)+1);
H.vox_offset  = typecast(blob(109:112),'single');
H.scl_scope   = typecast(blob(113:116),'single');
H.scl_inter   = typecast(blob(117:120),'single');
H.slice_end   = typecast(blob(121:122),'int16');
H.slice_code  = blob(123);
H.xyzt_units  = blob(124);
H.cal_max     = typecast(blob(125:128),'single');
H.cal_min     = typecast(blob(129:132),'single');
H.slice_duration = typecast(blob(133:136),'single');
H.toffset     = typecast(blob(137:140),'single');
H.glmax       = typecast(blob(141:144),'int32');
H.glmin       = typecast(blob(145:148),'int32');
H.descrip     = cstr2matlab(blob(149:228));
H.aux_file    = cstr2matlab(blob(229:252));
H.qform_code  = typecast(blob(253:254),'int16');
H.sform_code  = typecast(blob(255:256),'int16');
quats = typecast(blob(257:280),'single');
H.quatern_b = quats(1);
H.quatern_c = quats(2);
H.quatern_d = quats(3);
H.quatern_x = quats(4);
H.quatern_y = quats(5);
H.quatern_z = quats(6);
trafo = typecast(blob(281:328),'single');
H.srow_x = trafo(1:4);
H.srow_y = trafo(5:8);
H.srow_z = trafo(9:12);
%H.S = [H.srow_x; H.srow_y; H.srow_z; 0 0 0 1];
H.intent_name = cstr2matlab(blob(329:344));
H.magic = magic;


function ms = cstr2matlab(cs)
if cs(1)==0
  ms = '';
else
  ind = find(cs==0);
  if isempty(ind)
    ms = char(cs)';
  else 
    ms = char(cs(1:ind(1)-1))';
  end
end

   