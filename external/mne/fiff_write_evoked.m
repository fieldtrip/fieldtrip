function fiff_write_evoked(name,data)
%
% function fiff_write_evoked(name,data)
%
% name     filename
% data     the data structure returned from fiff_read_evoked
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
% Revision 1.14  2009/03/31 01:12:30  msh
% Improved ID handling
%
% Revision 1.13  2009/03/30 11:37:37  msh
% Added copying of measurement info blocks from the original like in mne_browse_raw
%
% Revision 1.12  2008/03/13 19:18:07  msh
% Read and write FIFF_MEAS_DATE from/to FIFFB_MEAS_INFO as appropriate
%
% Revision 1.11  2006/10/04 20:12:37  msh
% Added fiff_read_evoked_all
% Modified fiff_pick_channels_evoked to handle multiple data sets
% Added bad channel writing to fiff_write_evoked
%
% Revision 1.10  2006/09/26 23:20:23  msh
% Moved block ID and parent block id to the meas block.
%
% Revision 1.9  2006/06/14 05:32:43  msh
% Added missing FIFF_NAVE to each aspect written.
%
% Revision 1.8  2006/05/29 16:57:32  msh
% Make sure that scan numbers are correct.
%
% Revision 1.7  2006/04/23 15:29:40  msh
% Added MGH to the copyright
%
% Revision 1.6  2006/04/18 20:44:46  msh
% Added reading of forward solution.
% Use length instead of size when appropriate
%
% Revision 1.5  2006/04/14 15:49:49  msh
% Improved the channel selection code and added ch_names to measurement info.
%
% Revision 1.4  2006/04/14 11:03:57  msh
% Changed fiff_write_id write a given id.
% Added parent id writing.
%
% Revision 1.3  2006/04/13 23:09:46  msh
% Further streamlining of the coordinate transformations.
%
% Revision 1.2  2006/04/12 13:13:51  msh
% Added fiff_find_evoked.m
% Use aspect_kind field name instead of aspect_type
%
% Revision 1.1  2006/04/12 10:29:02  msh
% Made evoked data writing compatible with the structures returned in reading.
%
%
me='MNE:fiff_write_evoked';
if nargin ~= 2
    error(me,'File name and data required as an arguments');
end
%
global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end
%
%  Create the file and save the essentials
%
fid = fiff_start_file(name);
fiff_start_block(fid,FIFF.FIFFB_MEAS);
fiff_write_id(fid,FIFF.FIFF_BLOCK_ID);
if ~isempty(data.info.meas_id)
    fiff_write_id(fid,FIFF.FIFF_PARENT_BLOCK_ID,data.info.meas_id);
end
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
if length(blocks) > 0 && isfield(data.info,'filename') && ~isempty(data.info.filename)
    [ fid2, tree ] = fiff_open(data.info.filename);
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
%    General
%
fiff_write_float(fid,FIFF.FIFF_SFREQ,data.info.sfreq);
fiff_write_float(fid,FIFF.FIFF_HIGHPASS,data.info.highpass);
fiff_write_float(fid,FIFF.FIFF_LOWPASS,data.info.lowpass);
fiff_write_int(fid,FIFF.FIFF_NCHAN,data.info.nchan);
if [ ~isempty(data.info.meas_date) ]
    fiff_write_int(fid,FIFF.FIFF_MEAS_DATE,data.info.meas_date);
end
%
%    Coordinate transformations if the HPI result block was not there
%
if ~have_hpi_result
    if ~isempty(data.info.dev_head_t)
        fiff_write_coord_trans(fid,data.info.dev_head_t);
    end
    if ~isempty(data.info.ctf_head_t)
        fiff_write_coord_trans(fid,data.info.ctf_head_t);
    end
end
%
%  Channel information
%
for k = 1:data.info.nchan
    %
    %   Scan numbers may have been messed up
    %
    data.info.chs(k).scanno = k;
    fiff_write_ch_info(fid,data.info.chs(k));
end
%
%    Polhemus data
%
if ~isempty(data.info.dig) && ~have_isotrak
    fiff_start_block(fid,FIFF.FIFFB_ISOTRAK);
    for k = 1:length(data.info.dig)
        fiff_write_dig_point(fid,data.info.dig(k))
    end
    fiff_end_block(fid,FIFF.FIFFB_ISOTRAK);
end
%
%    Projectors
%
fiff_write_proj(fid,data.info.projs);
%
%    CTF compensation info
%
fiff_write_ctf_comp(fid,data.info.comps);
%
%    Bad channels
%
if length(data.info.bads) > 0
    fiff_start_block(fid,FIFF.FIFFB_MNE_BAD_CHANNELS);
    fiff_write_name_list(fid,FIFF.FIFF_MNE_CH_NAME_LIST,data.info.bads);
    fiff_end_block(fid,FIFF.FIFFB_MNE_BAD_CHANNELS);
end
%
%
fiff_end_block(fid,FIFF.FIFFB_MEAS_INFO);
%
% One or more evoked data sets
%
fiff_start_block(fid,FIFF.FIFFB_PROCESSED_DATA);
for set = 1:length(data.evoked)
    fiff_start_block(fid,FIFF.FIFFB_EVOKED);
    %
    % Comment is optional
    %
    if size(data.evoked(set).comment,2) > 0
        fiff_write_string(fid,FIFF.FIFF_COMMENT,data.evoked(set).comment);
    end
    %
    % First and last sample
    %
    fiff_write_int(fid,FIFF.FIFF_FIRST_SAMPLE,data.evoked(set).first);
    fiff_write_int(fid,FIFF.FIFF_LAST_SAMPLE,data.evoked(set).last);
    %
    % The epoch itself
    %
    fiff_start_block(fid,FIFF.FIFFB_ASPECT);
    %
    fiff_write_int(fid,FIFF.FIFF_ASPECT_KIND, ...
        data.evoked(set).aspect_kind);
    fiff_write_int(fid,FIFF.FIFF_NAVE,data.evoked(set).nave);
    decal = zeros(data.info.nchan,data.info.nchan);
    for k = 1:data.info.nchan
        decal(k,k) = 1.0/(data.info.chs(k).cal);
    end
    fiff_write_float_matrix(fid,FIFF.FIFF_EPOCH,decal*data.evoked(set).epochs);
    %
    fiff_end_block(fid,FIFF.FIFFB_ASPECT);
    %
    fiff_end_block(fid,FIFF.FIFFB_EVOKED);
end
fiff_end_block(fid,FIFF.FIFFB_PROCESSED_DATA);

fiff_end_block(fid,FIFF.FIFFB_MEAS);

fiff_end_file(fid);
return;
