function [fid,cals] = fiff_start_writing_raw(name,info,sel)
%
% function [fid,cals] = fiff_start_writing_raw(name,info,sel)
%
% name       filename
% info       The measurement info block of the source file
% sel        Which channels will be included in the output file (optional)
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
% Revision 1.6  2009/03/31 01:12:30  msh
% Improved ID handling
%
% Revision 1.5  2009/03/30 11:37:37  msh
% Added copying of measurement info blocks from the original like in mne_browse_raw
%
% Revision 1.4  2008/05/06 20:40:56  msh
% Fixed ordering of output for compatibility with maxfilter averager
%
% Revision 1.3  2008/04/16 22:24:57  msh
% Added megacq parameters to the measurement info
%
% Revision 1.2  2008/03/13 19:18:06  msh
% Read and write FIFF_MEAS_DATE from/to FIFFB_MEAS_INFO as appropriate
%
% Revision 1.1  2007/11/07 16:05:05  msh
% New routines for writing raw files
%

me='MNE:fiff_start_writing_raw';
if nargin ~= 2 && nargin ~= 3
    error(me,'Incorrect number of arguments');
end
%
global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end
%
%   We will always write floats
%
data_type = 4;
if nargin == 2
    sel = 1:info.nchan;
end
chs = info.chs(sel);
nchan = length(chs);
%
%  Create the file and save the essentials
%
fid = fiff_start_file(name);
fiff_start_block(fid,FIFF.FIFFB_MEAS);
fiff_write_id(fid,FIFF.FIFF_BLOCK_ID);
if ~isempty(info.meas_id)
    fiff_write_id(fid,FIFF.FIFF_PARENT_BLOCK_ID,info.meas_id);
end
%
%
%    Measurement info
%
fiff_start_block(fid,FIFF.FIFFB_MEAS_INFO);
%
%    Blocks from the original
%
blocks = [ FIFF.FIFFB_SUBJECT FIFF.FIFFB_HPI_MEAS FIFF.FIFFB_HPI_RESULT FIFF.FIFFB_ISOTRAK FIFF.FIFFB_PROCESSING_HISTORY ];
have_hpi_result = false;
have_isotrak    = false;
if length(blocks) > 0 && isfield(info,'filename') && ~isempty(info.filename)
    [ fid2, tree ] = fiff_open(info.filename);
    for k = 1:length(blocks)
        nodes = fiff_dir_tree_find(tree,blocks(k));
        fiff_copy_tree(fid2,tree.id,nodes,fid);
        if blocks(k) == FIFF.FIFFB_HPI_RESULT && length(nodes) > 0
            have_hpi_result = true;
        end
        if blocks(k) == FIFF.FIFFB_ISOTRAK && length(nodes) > 0
            have_isotrak = true;
        end
    end
    fclose(fid2);
end
%
%    megacq parameters
%
if ~isempty(info.acq_pars) || ~isempty(info.acq_stim)
    fiff_start_block(fid,FIFF.FIFFB_DACQ_PARS);
    if ~isempty(info.acq_pars)
        fiff_write_string(fid,FIFF.FIFF_DACQ_PARS, ...
            info.acq_pars);
    end
    if ~isempty(info.acq_stim)
        fiff_write_string(fid,FIFF.FIFF_DACQ_STIM, ...
            info.acq_stim);
    end
    fiff_end_block(fid,FIFF.FIFFB_DACQ_PARS);
end
%
%    Coordinate transformations if the HPI result block was not there
%
if ~have_hpi_result
    if ~isempty(info.dev_head_t)
        fiff_write_coord_trans(fid,info.dev_head_t);
    end
    if ~isempty(info.ctf_head_t)
        fiff_write_coord_trans(fid,info.ctf_head_t);
    end
end
%
%    Polhemus data
%
if ~isempty(info.dig) && ~have_isotrak
    fiff_start_block(fid,FIFF.FIFFB_ISOTRAK);
    for k = 1:length(info.dig)
        fiff_write_dig_point(fid,info.dig(k))
    end
    fiff_end_block(fid,FIFF.FIFFB_ISOTRAK);
end
%
%    Projectors
%
fiff_write_proj(fid,info.projs);
%
%    CTF compensation info
%
fiff_write_ctf_comp(fid,info.comps);
%
%    Bad channels
%
if length(info.bads) > 0
    fiff_start_block(fid,FIFF.FIFFB_MNE_BAD_CHANNELS);
    fiff_write_name_list(fid,FIFF.FIFF_MNE_CH_NAME_LIST,info.bads);
    fiff_end_block(fid,FIFF.FIFFB_MNE_BAD_CHANNELS);
end
%
%    General
%
fiff_write_float(fid,FIFF.FIFF_SFREQ,info.sfreq);
fiff_write_float(fid,FIFF.FIFF_HIGHPASS,info.highpass);
fiff_write_float(fid,FIFF.FIFF_LOWPASS,info.lowpass);
fiff_write_int(fid,FIFF.FIFF_NCHAN,nchan);
fiff_write_int(fid,FIFF.FIFF_DATA_PACK,data_type);
if [ ~isempty(info.meas_date) ]
    fiff_write_int(fid,FIFF.FIFF_MEAS_DATE,info.meas_date);
end
%
%    Channel info
%
for k = 1:nchan
    %
    %   Scan numbers may have been messed up
    %
    chs(k).scanno = k;
    chs(k).range  = 1.0;
    cals(k) = chs(k).cal;
    fiff_write_ch_info(fid,chs(k));
end
%
%
fiff_end_block(fid,FIFF.FIFFB_MEAS_INFO);
%
% Start the raw data
%
fiff_start_block(fid,FIFF.FIFFB_RAW_DATA);

return;
