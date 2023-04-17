function fiff_write_raw_buffer(fid,buf,cals,datatype)
%
% function fiff_write_raw_buffer(fid,buf,cals,datatype)
%
% fid        of an open raw data file
% buf        the buffer to write
% cals       calibration factors
% datatype   (optional) datatype to write, default float
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
% Revision 1.1  2007/11/07 16:05:05  msh
% New routines for writing raw files
%

me='MNE:fiff_write_raw_buffer';
if nargin < 3
    error(me,'Incorrect number of arguments');
end

if size(buf,1) ~= length(cals)
    error(me,'buffer and calibration sizes do not match');
end

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin < 4
    datatype = FIFF.FIFFT_FLOAT;
end

if datatype ~= FIFF.FIFFT_FLOAT
    warning(me, 'reading and writing of data in numeric precision ~= float is only supported in FieldTrip and MNE-Python');
end

switch datatype
    case FIFF.FIFFT_FLOAT
        fiff_write_float(fid,FIFF.FIFF_DATA_BUFFER,diag(1./cals)*buf);
    case FIFF.FIFFT_DOUBLE
        fiff_write_double(fid,FIFF.FIFF_DATA_BUFFER,diag(1./cals)*buf);
    case FIFF.FIFFT_COMPLEX_FLOAT
        fiff_write_complex(fid,FIFF.FIFF_DATA_BUFFER,diag(1./cals)*buf);
    case FIFF.FIFFT_COMPLEX_DOUBLE
        fiff_write_double_complex(fid,FIFF.FIFF_DATA_BUFFER,diag(1./cals)*buf);
    otherwise
        error(me,'unsupported datatype requested for writing of the buffer');
end

return;
