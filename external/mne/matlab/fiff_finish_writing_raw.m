function fiff_finish_writing_raw(fid)
%
% function fiff_finish_writing_raw(fid)
%
% fid        of an open raw data file
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
% Revision 1.2  2008/05/08 13:35:09  msh
% Ending FIFFB_MEAS was missing
%
% Revision 1.1  2007/11/07 16:05:05  msh
% New routines for writing raw files
%

me = 'MNE:fiff_finish_writing_raw';
if nargin ~= 1
    error(me, 'File id required as an argument');
end

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end
fiff_end_block(fid, FIFF.FIFFB_RAW_DATA);
fiff_end_block(fid, FIFF.FIFFB_MEAS);
fiff_end_file(fid);

return
