function blob = encode_nifti1(H)
%function blob = encode_nifti1(H)
%
% Encodes a NIFTI-1 header (=> raw 348 bytes (uint8)) from a Matlab structure
% that matches the C struct defined in nifti1.h.
%
% WARNING: This function currently ignores endianness !!!

% (C) 2010 S.Klanke

blob = uint8(zeros(1,348));

if ~isstruct(H)
	error 'Input must be a structure';
end
% see nift1.h for information on structure

sizeof_hdr = int32(348);
blob(1:4) = typecast(sizeof_hdr, 'uint8');

blob = setString(blob,  5, 14, H, 'data_type');
blob = setString(blob, 15, 32, H, 'db_name');
blob = setInt32( blob, 33, 36, H, 'extents');
blob = setInt16( blob, 37, 38, H, 'session_error');
blob = setInt8(  blob, 39, 39, H, 'regular');
blob = setInt8(  blob, 40, 40, H, 'dim_info');

dim = int16(H.dim(:)');
ndim = numel(dim);
if ndim<1 || ndim>7
	error 'Field "dim" must have 1..7 elements';
end
dim = [int16(ndim) dim];
blob(41:(42+2*ndim)) = typecast(dim,'uint8');

blob = setSingle(blob, 57, 60, H, 'intent_p1');
blob = setSingle(blob, 61, 64, H, 'intent_p2');
blob = setSingle(blob, 65, 68, H, 'intent_p3');
blob = setInt16( blob, 69, 70, H, 'intent_code');
blob = setInt16( blob, 71, 72, H, 'datatype');
blob = setInt16( blob, 73, 74, H, 'bitpix');
blob = setInt16( blob, 75, 76, H, 'slice_start');

blob = setSingle(blob, 77, 80, H, 'qfac');

if isfield(H,'pixdim')
	pixdim = single(H.pixdim(:)');
	ndim = numel(pixdim);
	if ndim<1 || ndim>7
		error 'Field "pixdim" must have 1..7 elements';
	end
	blob(81:(80+4*ndim)) = typecast(pixdim,'uint8');
end

blob = setSingle(blob, 109, 112, H, 'vox_offset');
blob = setSingle(blob, 113, 116, H, 'scl_scope');
blob = setSingle(blob, 117, 120, H, 'scl_inter');
blob = setInt16( blob, 121, 122, H, 'slice_end');
blob = setInt8(  blob, 123, 123, H, 'slice_code');
blob = setInt8(  blob, 124, 124, H, 'xyzt_units');
blob = setSingle(blob, 125, 128, H, 'cal_max');
blob = setSingle(blob, 129, 132, H, 'cal_min');
blob = setSingle(blob, 133, 136, H, 'slice_duration');
blob = setSingle(blob, 137, 140, H, 'toffset');
blob = setInt32( blob, 141, 144, H, 'glmax');
blob = setInt32( blob, 145, 148, H, 'glmin');
blob = setString(blob, 149, 228, H, 'descrip');
blob = setString(blob, 229, 252, H, 'aux_file');

blob = setInt16( blob, 253, 254, H, 'qform_code');
blob = setInt16( blob, 255, 256, H, 'sform_code');

blob = setSingle(blob, 257, 260, H, 'quatern_b');
blob = setSingle(blob, 261, 264, H, 'quatern_c');
blob = setSingle(blob, 265, 268, H, 'quatern_d');
blob = setSingle(blob, 269, 272, H, 'quatern_x');
blob = setSingle(blob, 273, 276, H, 'quatern_y');
blob = setSingle(blob, 277, 280, H, 'quatern_z');
blob = setSingle(blob, 281, 296, H, 'srow_x');
blob = setSingle(blob, 297, 312, H, 'srow_y');
blob = setSingle(blob, 313, 328, H, 'srow_z');

blob = setString(blob, 329, 344, H, 'intent_name');
if ~isfield(H,'magic')
	blob(345:347) = uint8('ni1');
else
	blob = setString(blob, 345, 347, H, 'magic');
end


function blob = setString(blob, begidx, endidx, H, fieldname)

if ~isfield(H,fieldname)
	return
end

F = getfield(H, fieldname);
ne = numel(F);
mx = endidx - begidx +1;
if ne > 0
	if ~ischar(F) || ne > mx
		errmsg = sprintf('Field "data_type" must be a string of maximally %i characters.', mx);
		error(errmsg);
	end
	blob(begidx:(begidx+ne-1)) = uint8(F(:)');
end


% set 32-bit integers (check #elements)
function blob = setInt32(blob, begidx, endidx, H, fieldname)
if ~isfield(H,fieldname)
	return
end

F = int32(getfield(H, fieldname));
ne = numel(F);
sp = (endidx - begidx +1) / 4;
if ne~=sp
	errmsg = sprintf('Field "data_type" must be an array with exactly %i elements.', sp);
	error(errmsg);
end
blob(begidx:(begidx+4*ne-1)) = typecast(F(:)', 'uint8');


% set 16-bit integers (check #elements)
function blob = setInt16(blob, begidx, endidx, H, fieldname)
if ~isfield(H,fieldname)
	return
end

F = int16(getfield(H, fieldname));
ne = numel(F);
sp = (endidx - begidx +1) / 2;
if ne~=sp
	errmsg = sprintf('Field "data_type" must be an array with exactly %i elements.', sp);
	error(errmsg);
end
blob(begidx:(begidx+2*ne-1)) = typecast(F(:)', 'uint8');


% just 8-bit integers (check #elements)
function blob = setInt8(blob, begidx, endidx, H, fieldname)
if ~isfield(H,fieldname)
	return
end

F = int8(getfield(H, fieldname));
ne = numel(F);
sp = (endidx - begidx +1);
if ne~=sp
	errmsg = sprintf('Field "data_type" must be an array with exactly %i elements.', sp);
	error(errmsg);
end
blob(begidx:(begidx+ne-1)) = typecast(F(:)', 'uint8');


% single precision floats
function blob = setSingle(blob, begidx, endidx, H, fieldname)
if ~isfield(H,fieldname)
	return
end

F = single(getfield(H, fieldname));
ne = numel(F);
sp = (endidx - begidx +1) / 4;
if ne~=sp
	errmsg = sprintf('Field "data_type" must be an array with exactly %i elements.', sp);
	error(errmsg);
end
blob(begidx:(begidx+4*ne-1)) = typecast(F(:)', 'uint8');



