function fiff_write_proj(fid,projs)
%
% fiff_write_proj(fid,projs)
%
% Writes the projection data into a fif file
%
%     fid           An open fif file descriptor
%     projs         The compensation data to write
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.9  2008/08/06 17:31:11  msh
%   Removed some debug output
%
%   Revision 1.8  2008/08/06 17:29:55  msh
%   Fixed missing end in writing the projection time.
%
%   Revision 1.7  2008/05/09 11:02:09  msh
%   Added FIFF_PROJ_ITEM_TIME for Neuromag compatibility
%
%   Revision 1.6  2008/05/06 20:40:56  msh
%   Fixed ordering of output for compatibility with maxfilter averager
%
%   Revision 1.5  2006/05/03 18:53:05  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.4  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.3  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.2  2006/04/12 10:51:18  msh
%   Added projection writing and compensation routines
%
%   Revision 1.1  2006/04/12 10:29:03  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%

me='MNE:fiff_write_proj';

if nargin ~= 2
    error(me,'Incorrect number of arguments');
end

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if isempty(projs)
    return;
end

fiff_start_block(fid,FIFF.FIFFB_PROJ);
for k = 1:length(projs)
    fiff_start_block(fid,FIFF.FIFFB_PROJ_ITEM);
    fiff_write_string(fid,FIFF.FIFF_NAME,projs(k).desc);
    fiff_write_int(fid,FIFF.FIFF_PROJ_ITEM_KIND,projs(k).kind);
    if projs(k).kind == FIFF.FIFFV_PROJ_ITEM_FIELD
        fiff_write_float(fid,FIFF.FIFF_PROJ_ITEM_TIME,0.0);
    end
    fiff_write_int(fid,FIFF.FIFF_NCHAN,projs(k).data.ncol);
    fiff_write_int(fid,FIFF.FIFF_PROJ_ITEM_NVEC, ...
        projs(k).data.nrow);
    fiff_write_int(fid,FIFF.FIFF_MNE_PROJ_ITEM_ACTIVE, ...
        projs(k).active);
    fiff_write_name_list(fid, ...
        FIFF.FIFF_PROJ_ITEM_CH_NAME_LIST, ...
        projs(k).data.col_names);
    fiff_write_float_matrix(fid,FIFF.FIFF_PROJ_ITEM_VECTORS,projs(k).data.data);
    fiff_end_block(fid,FIFF.FIFFB_PROJ_ITEM);
end
fiff_end_block(fid,FIFF.FIFFB_PROJ);

return;


end
