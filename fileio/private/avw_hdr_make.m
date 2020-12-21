function [ avw ] = avw_hdr_make

% AVW_HDR_MAKE - Create Analyze format data header (avw.hdr)
% 
% [ avw ] = avw_hdr_make
% 
% avw.hdr - a struct, all fields returned from the header.
%           For details, find a good description on the web
%           or see the Analyze File Format pdf in the 
%           mri_toolbox doc folder or see avw_hdr_read.m
% 
% See also, AVW_HDR_READ AVW_HDR_WRITE
%           AVW_IMG_READ AVW_IMG_WRITE
%

% $Revision$ $Date: 2009/01/14 09:24:45 $

% Licence:  GNU GPL, no express or implied warranties
% History:  06/2002, Darren.Weber@flinders.edu.au
%           02/2003, Darren.Weber@flinders.edu.au
%                    date/time bug at lines 97-98
%                    identified by Bennett.Landman@ieee.org
%
%                    The Analyze format is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '????????????????';
fprintf('\nAVW_HDR_MAKE [v%s]\n',version(12:16));  tic;


% Comments
% The header format is flexible and can be extended for new 
% user-defined data types. The essential structures of the header 
% are the header_key and the image_dimension.  See avw_hdr_read
% for more detail of the header structure
avw.hdr = make_header;

t=toc; fprintf('...done (%5.2f sec).\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ hdr ] = make_header
    
    hdr.hk   = header_key;
    hdr.dime = image_dimension;
    hdr.hist = data_history;
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hk] = header_key
    
    hk.sizeof_hdr       = int32(348);   % must be 348!
    hk.data_type(1:10)  = sprintf('%10s','');
    hk.db_name(1:18)    = sprintf('%18s','');
    hk.extents          = int32(16384);
    hk.session_error    = int16(0);
    hk.regular          = sprintf('%1s','r');  % might be uint8
    hk.hkey_un0         = uint8(0);
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dime ] = image_dimension
    
    dime.dim(1:8)       = int16([4 256 256 256 1 0 0 0]);
    dime.vox_units(1:4) = sprintf('%4s','mm');
    dime.cal_units(1:8) = sprintf('%8s','');
    dime.unused1        = int16(0);
    dime.datatype       = int16(2);
    dime.bitpix         = int16(8);
    dime.dim_un0        = int16(0);
    dime.pixdim(1:8)    = single([0 1 1 1 1000 0 0 0]);
    dime.vox_offset     = single(0);
    dime.funused1       = single(0);
    dime.funused2       = single(0);
    % Set default 8bit intensity scale (from MRIcro), otherwise funused3
    dime.roi_scale      = single(1);
    dime.cal_max        = single(0);
    dime.cal_min        = single(0);
    dime.compressed     = int32(0);
    dime.verified       = int32(0);
    dime.glmax          = int32(255);
    dime.glmin          = int32(0);
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ hist ] = data_history
    
    datime = clock;
    
    hist.descrip(1:80)      = sprintf('%-80s','mri_toolbox @ http://eeg.sf.net/');
    hist.aux_file(1:24)     = sprintf('%-24s','');
    hist.orient             = uint8(0); % sprintf( '%1s',''); % see notes in avw_hdr_read
    hist.originator(1:10)   = sprintf('%-10s','');
    hist.generated(1:10)    = sprintf('%-10s','mri_toolbx');
    hist.scannum(1:10)      = sprintf('%-10s','');
    hist.patient_id(1:10)   = sprintf('%-10s','');
    hist.exp_date(1:10)     = sprintf('%02d-%02d-%04d',datime(3),datime(2),datime(1));
    hist.exp_time(1:10)     = sprintf('%02d-%02d-%04.1f',datime(4),datime(5),datime(6));
    hist.hist_un0(1:3)      = sprintf( '%-3s','');
    hist.views              = int32(0);
    hist.vols_added         = int32(0);
    hist.start_field        = int32(0);
    hist.field_skip         = int32(0);
    hist.omax               = int32(0);
    hist.omin               = int32(0);
    hist.smax               = int32(0);
    hist.smin               = int32(0);
    
return
