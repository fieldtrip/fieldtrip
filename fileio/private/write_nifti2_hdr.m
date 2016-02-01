function write_nifti2_hdr(filename, hdr)

% WRITE_NIFTI2_HDR
%
% Use as
%   write_nifti2_hdr(filename, hdr)
% where
%   filename   = string
%   hdr        = structure with nifti-2 header information
%
% See also READ_NIFTI_HDR, READ_CIFTI, WRITE_CIFTI

% Copyright (C) 2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if ischar(filename)
  fid = fopen(filename, 'w');
else
  fid = filename;
end

fwrite(fid, 540, 'int32');

assert(fwrite(fid, hdr.magic          , 'int8'   )==8 ); % 4       `n', '+', `2', `\0','\r','\n','\032','\n' or (0x6E,0x2B,0x32,0x00,0x0D,0x0A,0x1A,0x0A)
assert(fwrite(fid, hdr.datatype       , 'int16'  )==1 ); % 12      See file formats
assert(fwrite(fid, hdr.bitpix         , 'int16'  )==1 ); % 14      See file formats
assert(fwrite(fid, hdr.dim            , 'int64'  )==8 ); % 16      See file formats
assert(fwrite(fid, hdr.intent_p1      , 'double' )==1 ); % 80      0
assert(fwrite(fid, hdr.intent_p2      , 'double' )==1 ); % 88      0
assert(fwrite(fid, hdr.intent_p3      , 'double' )==1 ); % 96      0
assert(fwrite(fid, hdr.pixdim         , 'double' )==8 ); % 104     0,1,1,1,1,1,1,1
assert(fwrite(fid, hdr.vox_offset     , 'int64'  )==1 ); % 168     Offset of data, minimum=544
assert(fwrite(fid, hdr.scl_slope      , 'double' )==1 ); % 176     1
assert(fwrite(fid, hdr.scl_inter      , 'double' )==1 ); % 184     0
assert(fwrite(fid, hdr.cal_max        , 'double' )==1 ); % 192     0
assert(fwrite(fid, hdr.cal_min        , 'double' )==1 ); % 200     0
assert(fwrite(fid, hdr.slice_duration , 'double' )==1 ); % 208     0
assert(fwrite(fid, hdr.toffset        , 'double' )==1 ); % 216     0
assert(fwrite(fid, hdr.slice_start    , 'int64'  )==1 ); % 224     0
assert(fwrite(fid, hdr.slice_end      , 'int64'  )==1 ); % 232     0
assert(fwrite(fid, hdr.descrip        , 'int8'   )==80); % 240     All zeros
assert(fwrite(fid, hdr.aux_file       , 'int8'   )==24); % 320     All zeros
assert(fwrite(fid, hdr.qform_code     , 'int32'  )==1 ); % 344     NIFTI_XFORM_UNKNOWN (0)
assert(fwrite(fid, hdr.sform_code     , 'int32'  )==1 ); % 348     NIFTI_XFORM_UNKNOWN (0)
assert(fwrite(fid, hdr.quatern_b      , 'double' )==1 ); % 352     0
assert(fwrite(fid, hdr.quatern_c      , 'double' )==1 ); % 360     0
assert(fwrite(fid, hdr.quatern_d      , 'double' )==1 ); % 368     0
assert(fwrite(fid, hdr.qoffset_x      , 'double' )==1 ); % 376     0
assert(fwrite(fid, hdr.qoffset_y      , 'double' )==1 ); % 384     0
assert(fwrite(fid, hdr.qOffset_z      , 'double' )==1 ); % 392     0
assert(fwrite(fid, hdr.srow_x         , 'double' )==4 ); % 400     0,0,0,0
assert(fwrite(fid, hdr.srow_y         , 'double' )==4 ); % 432     0,0,0,0
assert(fwrite(fid, hdr.srow_z         , 'double' )==4 ); % 464     0,0,0,0
assert(fwrite(fid, hdr.slice_code     , 'int32'  )==1 ); % 496     0
assert(fwrite(fid, hdr.xyzt_units     , 'int32'  )==1 ); % 500     0xC (seconds, millimeters)
assert(fwrite(fid, hdr.intent_code    , 'int32'  )==1 ); % 504     See file formats
assert(fwrite(fid, hdr.intent_name    , 'int8'   )==16); % 508     See file formats
assert(fwrite(fid, hdr.dim_info       , 'int8'   )==1 ); % 524     0
assert(fwrite(fid, hdr.unused_str     , 'int8'   )==15); % 525     All zeros
% disp(ftell(fid));                                     % 540     End of the header

if ischar(filename)
  fclose(fid);
end

