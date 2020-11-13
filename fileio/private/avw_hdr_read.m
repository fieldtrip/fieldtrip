function [ avw, machine ] = avw_hdr_read(fileprefix, machine, verbose)

% avw_hdr_read - read Analyze format data header (*.hdr)
%
% [ avw, machine ] = avw_hdr_read(fileprefix, [machine], [verbose])
%
% fileprefix - string filename (without .hdr); the file name
%              can be given as a full path or relative to the
%              current directory.
%
% machine - a string, see machineformat in fread for details.
%           The default here is 'ieee-le' but the routine
%           will automatically switch between little and big
%           endian to read any such Analyze header.  It
%           reports the appropriate machine format and can
%           return the machine value.
%
% avw.hdr - a struct, all fields returned from the header.
%           For details, find a good description on the web
%           or see the Analyze File Format pdf in the
%           mri_toolbox doc folder or read this .m file.
%
% verbose - the default is to output processing information to the command
%           window.  If verbose = 0, this will not happen.
%
% This function is called by avw_img_read
%
% See also avw_hdr_write, avw_hdr_make, avw_view_hdr, avw_view
%

% $Revision$ $Date: 2009/01/14 09:24:45 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%                    The Analyze format and c code below is copyright
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('verbose','var'), verbose = 1; end

if verbose,
    version = '[$Revision$]';
    fprintf('\nAVW_HDR_READ [v%s]\n',version(12:16));  tic;
end

if ~exist('fileprefix','var'),
    ft_error('...no input fileprefix - see help avw_hdr_read\n\n');
end
if ~exist('machine','var'), machine = 'ieee-le'; end


if contains(fileprefix, '.hdr')
    % fprintf('...removing .hdr extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.hdr','');
end
if contains(fileprefix, '.img')
    % fprintf('...removing .img extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.img','');
end
file = sprintf('%s.hdr',fileprefix);

if exist(file),
    if verbose,
        fprintf('...reading %s Analyze format',machine);
    end
    fid = fopen_or_error(file,'r',machine);
    avw.hdr = read_header(fid,verbose);
    avw.fileprefix = fileprefix;
    fclose(fid);

    if ~isequal(avw.hdr.hk.sizeof_hdr,348),
        if verbose, fprintf('...failed.\n'); end
        % first try reading the opposite endian to 'machine'
        switch machine,
        case 'ieee-le', machine = 'ieee-be';
        case 'ieee-be', machine = 'ieee-le';
        end
        if verbose, fprintf('...reading %s Analyze format',machine); end
        fid = fopen_or_error(file,'r',machine);
        avw.hdr = read_header(fid,verbose);
        avw.fileprefix = fileprefix;
        fclose(fid);
    end
    if ~isequal(avw.hdr.hk.sizeof_hdr,348)
        % Now throw an error
        if verbose, fprintf('...failed.\n'); end
        ft_error('...size of header not equal to 348 bytes!\n\n');
    end
else
    ft_error('...cannot find file %s.hdr\n\n',file);
end

if verbose
    t=toc; fprintf('...done (%5.2f sec).\n',t);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dsr ] = read_header(fid,verbose)

    % Original header structures - ANALYZE 7.5
    %struct dsr
    %       {
    %       struct header_key hk;            /*   0 +  40       */
    %       struct image_dimension dime;     /*  40 + 108       */
    %       struct data_history hist;        /* 148 + 200       */
    %       };                               /* total= 348 bytes*/
    dsr.hk   = header_key(fid);
    dsr.dime = image_dimension(fid,verbose);
    dsr.hist = data_history(fid);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hk] = header_key(fid)

    % The required elements in the header_key substructure are:
    %
    % int sizeof_header   Must indicate the byte size of the header file.
    % int extents         Should be 16384, the image file is created as
    %                     contiguous with a minimum extent size.
    % char regular        Must be 'r' to indicate that all images and
    %                     volumes are the same size.

    % Original header structures - ANALYZE 7.5
    % struct header_key                      /* header key      */
    %       {                                /* off + size      */
    %       int sizeof_hdr                   /*  0 +  4         */
    %       char data_type[10];              /*  4 + 10         */
    %       char db_name[18];                /* 14 + 18         */
    %       int extents;                     /* 32 +  4         */
    %       short int session_error;         /* 36 +  2         */
    %       char regular;                    /* 38 +  1         */
    %       char hkey_un0;                   /* 39 +  1         */
    %       };                               /* total=40 bytes  */

    fseek(fid,0,'bof');

    hk.sizeof_hdr    = fread(fid, 1,'*int32');  % should be 348!
    hk.data_type     = fread(fid,10,'*char')';
    hk.db_name       = fread(fid,18,'*char')';
    hk.extents       = fread(fid, 1,'*int32');
    hk.session_error = fread(fid, 1,'*int16');
    hk.regular       = fread(fid, 1,'*char')'; % might be uint8
    hk.hkey_un0      = fread(fid, 1,'*uint8')';

    % check if this value was a char zero
    if hk.hkey_un0 == 48,
        hk.hkey_un0 = 0;
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dime ] = image_dimension(fid,verbose)

    %struct image_dimension
    %       {                                /* off + size      */
    %       short int dim[8];                /* 0 + 16          */
    %           /*
    %           dim[0]      Number of dimensions in database; usually 4.
    %           dim[1]      Image X dimension;  number of *pixels* in an image row.
    %           dim[2]      Image Y dimension;  number of *pixel rows* in slice.
    %           dim[3]      Volume Z dimension; number of *slices* in a volume.
    %           dim[4]      Time points; number of volumes in database
    %           */
    %       char vox_units[4];               /* 16 + 4          */
    %       char cal_units[8];               /* 20 + 8          */
    %       short int unused1;               /* 28 + 2          */
    %       short int datatype;              /* 30 + 2          */
    %       short int bitpix;                /* 32 + 2          */
    %       short int dim_un0;               /* 34 + 2          */
    %       float pixdim[8];                 /* 36 + 32         */
    %           /*
    %               pixdim[] specifies the voxel dimensions:
    %               pixdim[1] - voxel width, mm
    %               pixdim[2] - voxel height, mm
    %               pixdim[3] - slice thickness, mm
    %               pixdim[4] - volume timing, in msec
    %                   ..etc
    %           */
    %       float vox_offset;                /* 68 + 4          */
    %       float roi_scale;                 /* 72 + 4          */
    %       float funused1;                  /* 76 + 4          */
    %       float funused2;                  /* 80 + 4          */
    %       float cal_max;                   /* 84 + 4          */
    %       float cal_min;                   /* 88 + 4          */
    %       int compressed;                  /* 92 + 4          */
    %       int verified;                    /* 96 + 4          */
    %       int glmax;                       /* 100 + 4         */
    %       int glmin;                       /* 104 + 4         */
    %       };                               /* total=108 bytes */

    dime.dim        = fread(fid,8,'*int16')';
    dime.vox_units  = fread(fid,4,'*char')';
    dime.cal_units  = fread(fid,8,'*char')';
    dime.unused1    = fread(fid,1,'*int16');
    dime.datatype   = fread(fid,1,'*int16');
    dime.bitpix     = fread(fid,1,'*int16');
    dime.dim_un0    = fread(fid,1,'*int16');
    dime.pixdim     = fread(fid,8,'*float')';
    dime.vox_offset = fread(fid,1,'*float');
    dime.roi_scale  = fread(fid,1,'*float');
    dime.funused1   = fread(fid,1,'*float');
    dime.funused2   = fread(fid,1,'*float');
    dime.cal_max    = fread(fid,1,'*float');
    dime.cal_min    = fread(fid,1,'*float');
    dime.compressed = fread(fid,1,'*int32');
    dime.verified   = fread(fid,1,'*int32');
    dime.glmax      = fread(fid,1,'*int32');
    dime.glmin      = fread(fid,1,'*int32');

    if dime.dim(1) < 4, % Number of dimensions in database; usually 4.
        if verbose,
            fprintf('...ensuring 4 dimensions in avw.hdr.dime.dim\n');
        end
        dime.dim(1) = int16(4);
    end
    if dime.dim(5) < 1, % Time points; number of volumes in database
        if verbose,
            fprintf('...ensuring at least 1 volume in avw.hdr.dime.dim(5)\n');
        end
        dime.dim(5) = int16(1);
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ hist ] = data_history(fid)

    % Original header structures - ANALYZE 7.5
    %struct data_history
    %       {                                /* off + size      */
    %       char descrip[80];                /* 0 + 80          */
    %       char aux_file[24];               /* 80 + 24         */
    %       char orient;                     /* 104 + 1         */
    %       char originator[10];             /* 105 + 10        */
    %       char generated[10];              /* 115 + 10        */
    %       char scannum[10];                /* 125 + 10        */
    %       char patient_id[10];             /* 135 + 10        */
    %       char exp_date[10];               /* 145 + 10        */
    %       char exp_time[10];               /* 155 + 10        */
    %       char hist_un0[3];                /* 165 + 3         */
    %       int views                        /* 168 + 4         */
    %       int vols_added;                  /* 172 + 4         */
    %       int start_field;                 /* 176 + 4         */
    %       int field_skip;                  /* 180 + 4         */
    %       int omax;                        /* 184 + 4         */
    %       int omin;                        /* 188 + 4         */
    %       int smax;                        /* 192 + 4         */
    %       int smin;                        /* 196 + 4         */
    %       };                               /* total=200 bytes */

    hist.descrip     = fread(fid,80,'*char')';
    hist.aux_file    = fread(fid,24,'*char')';
    hist.orient      = fread(fid, 1,'*uint8');  % see note below on char
    hist.originator  = fread(fid,10,'*char')';
    hist.generated   = fread(fid,10,'*char')';
    hist.scannum     = fread(fid,10,'*char')';
    hist.patient_id  = fread(fid,10,'*char')';
    hist.exp_date    = fread(fid,10,'*char')';
    hist.exp_time    = fread(fid,10,'*char')';
    hist.hist_un0    = fread(fid, 3,'*char')';
    hist.views       = fread(fid, 1,'*int32');
    hist.vols_added  = fread(fid, 1,'*int32');
    hist.start_field = fread(fid, 1,'*int32');
    hist.field_skip  = fread(fid, 1,'*int32');
    hist.omax        = fread(fid, 1,'*int32');
    hist.omin        = fread(fid, 1,'*int32');
    hist.smax        = fread(fid, 1,'*int32');
    hist.smin        = fread(fid, 1,'*int32');

    % check if hist.orient was saved as ascii char value
    switch hist.orient,
        case 48, hist.orient = uint8(0);
        case 49, hist.orient = uint8(1);
        case 50, hist.orient = uint8(2);
        case 51, hist.orient = uint8(3);
        case 52, hist.orient = uint8(4);
        case 53, hist.orient = uint8(5);
    end

return


% Note on using char:
% The 'char orient' field in the header is intended to
% hold simply an 8-bit unsigned integer value, not the ASCII representation
% of the character for that value.  A single 'char' byte is often used to
% represent an integer value in Analyze if the known value range doesn't
% go beyond 0-255 - saves a byte over a short int, which may not mean
% much in today's computing environments, but given that this format
% has been around since the early 1980's, saving bytes here and there on
% older systems was important!  In this case, 'char' simply provides the
% byte of storage - not an indicator of the format for what is stored in
% this byte.  Generally speaking, anytime a single 'char' is used, it is
% probably meant to hold an 8-bit integer value, whereas if this has
% been dimensioned as an array, then it is intended to hold an ASCII
% character string, even if that was only a single character.
% Denny  <hanson.dennis2@mayo.edu>


% Comments
% The header format is flexible and can be extended for new
% user-defined data types. The essential structures of the header
% are the header_key and the image_dimension.
%

% The required elements in the header_key substructure are:
%
% int sizeof_header   Must indicate the byte size of the header file.
% int extents         Should be 16384, the image file is created as
%                     contiguous with a minimum extent size.
% char regular        Must be 'r' to indicate that all images and
%                     volumes are the same size.
%

% The image_dimension substructure describes the organization and
% size of the images. These elements enable the database to reference
% images by volume and slice number. Explanation of each element follows:
%
% short int dim[ ];      /* Array of the image dimensions */
%
% dim[0]      Number of dimensions in database; usually 4.
% dim[1]      Image X dimension; number of pixels in an image row.
% dim[2]      Image Y dimension; number of pixel rows in slice.
% dim[3]      Volume Z dimension; number of slices in a volume.
% dim[4]      Time points; number of volumes in database.
% dim[5]      Undocumented.
% dim[6]      Undocumented.
% dim[7]      Undocumented.
%
% char vox_units[4]     Specifies the spatial units of measure for a voxel.
% char cal_units[8]      Specifies the name of the calibration unit.
% short int unused1      /* Unused */
% short int datatype      /* Datatype for this image set */
% /*Acceptable values for datatype are*/
% #define DT_NONE             0
% #define DT_UNKNOWN          0    /*Unknown data type*/
% #define DT_BINARY           1    /*Binary             ( 1 bit per voxel)*/
% #define DT_UNSIGNED_CHAR    2    /*Unsigned character ( 8 bits per voxel)*/
% #define DT_SIGNED_SHORT     4    /*Signed short       (16 bits per voxel)*/
% #define DT_SIGNED_INT       8    /*Signed integer     (32 bits per voxel)*/
% #define DT_FLOAT           16    /*Floating point     (32 bits per voxel)*/
% #define DT_COMPLEX         32    /*Complex (64 bits per voxel; 2 floating point numbers)/*
% #define DT_DOUBLE          64    /*Double precision   (64 bits per voxel)*/
% #define DT_RGB            128    /*A Red-Green-Blue datatype*/
% #define DT_ALL            255    /*Undocumented*/
%
% short int bitpix;    /* Number of bits per pixel; 1, 8, 16, 32, or 64. */
% short int dim_un0;   /* Unused */
%
% float pixdim[];     Parallel array to dim[], giving real world measurements in mm and ms.
%       pixdim[0];    Pixel dimensions?
%       pixdim[1];    Voxel width in mm.
%       pixdim[2];    Voxel height in mm.
%       pixdim[3];    Slice thickness in mm.
%       pixdim[4];    timeslice in ms (ie, TR in fMRI).
%       pixdim[5];    Undocumented.
%       pixdim[6];    Undocumented.
%       pixdim[7];    Undocumented.
%
% float vox_offset;   Byte offset in the .img file at which voxels start. This value can be
%                     negative to specify that the absolute value is applied for every image
%                     in the file.
%
% float roi_scale; Specifies the Region Of Interest scale?
% float funused1; Undocumented.
% float funused2; Undocumented.
%
% float cal_max; Specifies the upper bound of the range of calibration values.
% float cal_min; Specifies the lower bound of the range of calibration values.
%
% int compressed; Undocumented.
% int verified;   Undocumented.
%
% int glmax;    The maximum pixel value for the entire database.
% int glmin;    The minimum pixel value for the entire database.
%
%
