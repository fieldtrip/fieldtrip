function [hdr] = read_nifti2_hdr(filename)

% READ_NIFTI2_HDR
%
% Use as
%   [hdr] = read_nifti2_hdr(filename)
% where
%   filename   = string
%   
% This implements the format as described at
%   http://www.nitrc.org/forum/forum.php?thread_id=2148&forum_id=1941
%
% Please note that it is different from the suggested format described here
%   http://www.nitrc.org/forum/forum.php?thread_id=2070&forum_id=1941
% and
%   https://mail.nmr.mgh.harvard.edu/pipermail//freesurfer/2011-February/017482.html
% Notably, the unused fields have been removed and the size has been
% reduced from 560 to 540 bytes.
%
% See also WRITE_NIFTI_HDR, READ_CIFTI, WRITE_CIFTI

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

hdr.endian = 'l';
fid = fopen(filename, 'rb', hdr.endian);
hdr.sizeof_hdr = fread(fid, [1 1 ], 'int32=>int32'); % 0

if hdr.sizeof_hdr~=348 && hdr.sizeof_hdr~=540
  % try opening as big endian
  fclose(fid);
  hdr.endian = 'b';
  fid = fopen(filename, 'r', 'b');
  hdr.sizeof_hdr = fread(fid, [1 1 ], 'int32=>int32'); % 0
end

if hdr.sizeof_hdr~=348 && hdr.sizeof_hdr~=540
  fclose(fid);
  error('cannot open %s as nifti file, hdr size = %d, should be 348 or 540\n', filename, hdr.sizeof_hdr);
else
  % the file is now open with the appropriate little or big-endianness
end

if hdr.sizeof_hdr==348
  % the remainder of the code is for nifti-2 files
  error('%s seems to be a nifti-1 file', filename)
end

hdr.magic           = fread(fid, [1 8 ], 'int8=>int8'     ); % 4       `n', '+', `2', `\0','\r','\n','\032','\n' or (0x6E,0x2B,0x32,0x00,0x0D,0x0A,0x1A,0x0A)
hdr.datatype        = fread(fid, [1 1 ], 'int16=>int16'   ); % 12      See file formats
hdr.bitpix          = fread(fid, [1 1 ], 'int16=>int16'   ); % 14      See file formats
hdr.dim             = fread(fid, [1 8 ], 'int64=>double'  ); % 16      See file formats

if hdr.dim(1)<1 || hdr.dim(1)>7
  % see http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/dim.html
  error('inconsistent endianness in the header');
end

hdr.intent_p1       = fread(fid, [1 1 ], 'double=>double' ); % 80      0
hdr.intent_p2       = fread(fid, [1 1 ], 'double=>double' ); % 88      0
hdr.intent_p3       = fread(fid, [1 1 ], 'double=>double' ); % 96      0
hdr.pixdim          = fread(fid, [1 8 ], 'double=>double' ); % 104     0,1,1,1,1,1,1,1
hdr.vox_offset      = fread(fid, [1 1 ], 'int64=>int64'   ); % 168     Offset of data, minimum=544
hdr.scl_slope       = fread(fid, [1 1 ], 'double=>double' ); % 176     1
hdr.scl_inter       = fread(fid, [1 1 ], 'double=>double' ); % 184     0
hdr.cal_max         = fread(fid, [1 1 ], 'double=>double' ); % 192     0
hdr.cal_min         = fread(fid, [1 1 ], 'double=>double' ); % 200     0
hdr.slice_duration  = fread(fid, [1 1 ], 'double=>double' ); % 208     0
hdr.toffset         = fread(fid, [1 1 ], 'double=>double' ); % 216     0
hdr.slice_start     = fread(fid, [1 1 ], 'int64=>int64'   ); % 224     0
hdr.slice_end       = fread(fid, [1 1 ], 'int64=>int64'   ); % 232     0
hdr.descrip         = fread(fid, [1 80], 'int8=>char'     ); % 240     All zeros
hdr.aux_file        = fread(fid, [1 24], 'int8=>char'     ); % 320     All zeros
hdr.qform_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 344     NIFTI_XFORM_UNKNOWN (0)
hdr.sform_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 348     NIFTI_XFORM_UNKNOWN (0)
hdr.quatern_b       = fread(fid, [1 1 ], 'double=>double' ); % 352     0
hdr.quatern_c       = fread(fid, [1 1 ], 'double=>double' ); % 360     0
hdr.quatern_d       = fread(fid, [1 1 ], 'double=>double' ); % 368     0
hdr.qoffset_x       = fread(fid, [1 1 ], 'double=>double' ); % 376     0
hdr.qoffset_y       = fread(fid, [1 1 ], 'double=>double' ); % 384     0
hdr.qOffset_z       = fread(fid, [1 1 ], 'double=>double' ); % 392     0
hdr.srow_x          = fread(fid, [1 4 ], 'double=>double' ); % 400     0,0,0,0
hdr.srow_y          = fread(fid, [1 4 ], 'double=>double' ); % 432     0,0,0,0
hdr.srow_z          = fread(fid, [1 4 ], 'double=>double' ); % 464     0,0,0,0
hdr.slice_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 496     0
hdr.xyzt_units      = fread(fid, [1 1 ], 'int32=>int32'   ); % 500     0xC (seconds, millimeters)
hdr.intent_code     = fread(fid, [1 1 ], 'int32=>int32'   ); % 504     See file formats
hdr.intent_name     = fread(fid, [1 16], 'int8=>char'     ); % 508     See file formats
hdr.dim_info        = fread(fid, [1 1 ], 'int8=>int8'     ); % 524     0
hdr.unused_str      = fread(fid, [1 15], 'int8=>char'     ); % 525     All zeros
% disp(ftell(fid));                                          % 540     End of the header

fclose(fid);