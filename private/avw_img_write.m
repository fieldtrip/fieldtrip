function avw_img_write(avw, fileprefix, IMGorient, machine, verbose)

% avw_img_write - write Analyze image files (*.img)
% 
% avw_img_write(avw,fileprefix,[IMGorient],[machine],[verbose])
% 
% avw.img    - a 3D matrix of image data (double precision).
% avw.hdr    - a struct with image data parameters.  If
%              not empty, this function calls avw_hdr_write.
% 
% fileprefix - a string, the filename without the .img
%              extension. If empty, may use avw.fileprefix
% 
% IMGorient - optional int, force writing of specified 
%             orientation, with values:
% 
%   [],  if empty, will use avw.hdr.hist.orient field
%    0,  transverse/axial unflipped (default, radiological)
%    1,  coronal unflipped
%    2,  sagittal unflipped
%    3,  transverse/axial flipped, left to right
%    4,  coronal flipped, anterior to posterior
%    5,  sagittal flipped, superior to inferior
% 
% This function will set avw.hdr.hist.orient and write the 
% image data in a corresponding order.  This function is 
% in alpha development, so it has not been exhaustively 
% tested (07/2003). See avw_img_read for more information 
% and documentation on the orientation option.  
% Orientations 3-5 are NOT recommended!  They are part 
% of the Analyze format, but only used in Analyze
% for faster raster graphics during movies.
% 
% machine - a string, see machineformat in fread for details.
%           The default here is 'ieee-le'.
% 
% verbose - the default is to output processing information to the command
%           window.  If verbose = 0, this will not happen.
%
% Tip: to change the data type, set avw.hdr.dime.datatype to:
% 
%     1    Binary             (  1 bit  per voxel)
%     2    Unsigned character (  8 bits per voxel)
%     4    Signed short       ( 16 bits per voxel)
%     8    Signed integer     ( 32 bits per voxel)
%    16    Floating point     ( 32 bits per voxel)
%    32    Complex, 2 floats  ( 64 bits per voxel), not supported
%    64    Double precision   ( 64 bits per voxel)
%   128    Red-Green-Blue     (128 bits per voxel), not supported
% 
% See also: avw_write, avw_hdr_write, 
%           avw_read, avw_hdr_read, avw_img_read, avw_view
% 

% $Revision$ $Date: 2009/01/14 09:24:45 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%                    The Analyze format is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------
% Check inputs

if ~exist('avw','var')
    doc avw_img_write;
    error('...no input avw.');
elseif isempty(avw)
    error('...empty input avw.');
elseif ~isfield(avw,'img')
    error('...empty input avw.img');
end

if ~exist('fileprefix','var')
    if isfield(avw,'fileprefix')
        if ~isempty(avw.fileprefix)
            fileprefix = avw.fileprefix;
        else
            fileprefix = [];
        end
    else
        fileprefix = [];
    end
end
if isempty(fileprefix)
    [fileprefix, pathname, filterindex] = uiputfile('*.hdr','Specify an output Analyze .hdr file');
    if pathname, cd(pathname); end
    if ~fileprefix
        doc avw_img_write;
        error('no output .hdr file specified');
    end
end

if findstr('.hdr',fileprefix)
%    fprintf('AVW_IMG_WRITE: Removing .hdr extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.hdr','');
end
if findstr('.img',fileprefix)
%    fprintf('AVW_IMG_WRITE: Removing .img extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.img','');
end

if ~exist('IMGorient','var'), IMGorient = ''; end
if ~exist('machine','var'), machine = 'ieee-le'; end
if ~exist('verbose','var'), verbose = 1; end

if isempty(IMGorient), IMGorient = ''; end
if isempty(machine), machine = 'ieee-le'; end
if isempty(verbose), verbose = 1; end



%------------------------------------------------------------------------
% MAIN
if verbose
    version = '[$Revision$]';
    fprintf('\nAVW_IMG_WRITE [v%s]\n',version(12:16));  tic;
end

fid = fopen(sprintf('%s.img',fileprefix),'w',machine);
if fid < 0
    error('Cannot open file %s.img\n',fileprefix);
else
    avw = write_image(fid,avw,fileprefix,IMGorient,machine,verbose);
end

if verbose
    t=toc; fprintf('...done (%5.2f sec).\n\n',t);
end

% MUST write header after the image, to ensure any
% orientation changes during image write are saved
% in the header
avw_hdr_write(avw,fileprefix,machine,verbose);

return


function avw_hdr_write(avw, fileprefix, machine, verbose)

% AVW_HDR_WRITE - Write Analyze header file (*.hdr)
% 
% avw_hdr_write(avw,[fileprefix],[machine],[verbose])
% 
% eg, avw_hdr_write(avw,'test');
% 
% avw        - a struct with .hdr field, which itself is a struct,
%              containing all fields of an Analyze header.
%              For details, see avw_hdr_read.m
% 
% fileprefix - a string, the filename without the .hdr extension.
%              If empty, may use avw.fileprefix
% 
% machine    - a string, see machineformat in fread for details.
%              The default here is 'ieee-le'.
% 
% verbose - the default is to output processing information to the command
%           window.  If verbose = 0, this will not happen.
%
% See also, AVW_HDR_READ AVW_HDR_MAKE
%           AVW_IMG_READ AVW_IMG_WRITE
%

% $Revision$ $Date: 2009/01/14 09:24:45 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%           02/2003, Bennett.Landman@ieee.org
%                    - more specific data history var sizes
%                    - 02/2003 confirmed, Darren
%
%                    The Analyze format and c code below is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('verbose','var'), verbose = 1; end
if ~exist('machine','var'), machine = 'ieee-le'; end

if verbose,
    version = '[$Revision$]';
    fprintf('AVW_HDR_WRITE [v%s]\n',version(12:16));  tic;
end

%----------------------------------------------------------------------------
% Check inputs

if ~exist('avw','var')
    warning('...no input avw - calling avw_hdr_make\n');
    avw = avw_hdr_make;
elseif isempty(avw)
    warning('...empty input avw - calling avw_hdr_make\n');
    avw = avw_hdr_make;
elseif ~isfield(avw,'hdr')
    warning('...empty input avw.hdr - calling avw_hdr_make\n');
    avw = avw_hdr_make;
end
if ~isequal(avw.hdr.hk.sizeof_hdr,348)
    error('...avw.hdr.hk.sizeof_hdr must be 348!\n');
end

quit = 0;
if ~exist('fileprefix','var')
    if isfield(avw,'fileprefix')
        if ~isempty(avw.fileprefix)
            fileprefix = avw.fileprefix;
        else
            quit = 1;
        end
    else
        quit = 1;
    end
    if quit
        helpwin avw_hdr_write;
        error('...no input fileprefix - see help avw_hdr_write\n\n');
        return;
    end
end

if findstr('.hdr',fileprefix)
%    fprintf('AVW_HDR_WRITE: Removing .hdr extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.hdr','');
end



%----------------------------------------------------------------------------
% MAIN

if verbose, tic; end


% % force volume to 4D if necessary; conforms to AVW standard.  i lifted this
% % code from mri_toolbox/avw_hdr_read.m and modified it just a little bit
% % (using minDim = 4; etc)
% minDim = 4;
% currDim = double(avw.hdr.dime.dim(1));
% if ( currDim < minDim )
%    % fprintf( 'Warning %s: Forcing %d dimensions in avw.hdr.dime.dim\n', ...
%    %    mfilename, minDim );
%    avw.hdr.dime.dim(1) = int16(minDim);
%    avw.hdr.dime.dim(currDim+2:minDim+1) = int16(1);
%    avw.hdr.dime.pixdim(1) = int16(minDim);
%    avw.hdr.dime.pixdim(currDim+2:minDim+1) = int16(1);
% end;


fid = fopen(sprintf('%s.hdr',fileprefix),'w',machine);
if fid < 0
    error('Cannot write to file %s.hdr\n',fileprefix);
else
    if verbose, fprintf('...writing %s Analyze header.\n',machine); end
    write_header(fid,avw,verbose);
end

if verbose, t=toc; fprintf('...done (%5.2f sec).\n\n',t); end

return


%-----------------------------------------------------------------------------------

function avw = write_image(fid,avw,fileprefix,IMGorient,machine,verbose)

% short int bitpix;    /* Number of bits per pixel; 1, 8, 16, 32, or 64. */ 
% short int datatype   /* Datatype for this image set */
% /*Acceptable values for datatype are*/ 
% #define DT_NONE             0
% #define DT_UNKNOWN          0    /*Unknown data type*/ 
% #define DT_BINARY           1    /*Binary             ( 1 bit per voxel)*/ 
% #define DT_UNSIGNED_CHAR    2    /*Unsigned character ( 8 bits per voxel)*/ 
% #define DT_SIGNED_SHORT     4    /*Signed short       (16 bits per voxel)*/ 
% #define DT_SIGNED_INT       8    /*Signed integer     (32 bits per voxel)*/ 
% #define DT_FLOAT           16    /*Floating point     (32 bits per voxel)*/ 
% #define DT_COMPLEX         32    /*Complex,2 floats   (64 bits per voxel)/* 
% #define DT_DOUBLE          64    /*Double precision   (64 bits per voxel)*/ 
% #define DT_RGB            128    /*A Red-Green-Blue datatype*/
% #define DT_ALL            255    /*Undocumented*/

switch double(avw.hdr.dime.datatype)
case   1
    avw.hdr.dime.bitpix = int16( 1); precision = 'bit1';
case   2
    avw.hdr.dime.bitpix = int16( 8); precision = 'uchar';
case   4
    avw.hdr.dime.bitpix = int16(16); precision = 'int16';
case   8
    avw.hdr.dime.bitpix = int16(32); precision = 'int32';
case  16
    avw.hdr.dime.bitpix = int16(32); precision = 'single';
case  32
    error('...complex datatype not yet supported.\n');
case  64
    avw.hdr.dime.bitpix = int16(64); precision = 'double';
case 128
    error('...RGB datatype not yet supported.\n');
otherwise
    warning('...unknown datatype, using type 16 (32 bit floats).\n');
    avw.hdr.dime.datatype = int16(16);
    avw.hdr.dime.bitpix = int16(32); precision = 'single';
end


% write the .img file, depending on the .img orientation
if verbose,
    fprintf('...writing %s precision Analyze image (%s).\n',precision,machine);
end

fseek(fid,0,'bof');

% The standard image orientation is axial unflipped
if isempty(avw.hdr.hist.orient),
    msg = [ '...avw.hdr.hist.orient ~= 0.\n',...
            '   This function assumes the input avw.img is\n',...
            '   in axial unflipped orientation in memory.  This is\n',...
            '   created by the avw_img_read function, which converts\n',...
            '   any input file image to axial unflipped in memory.\n'];
    warning(msg)
end

if isempty(IMGorient),
    if verbose,
        fprintf('...no IMGorient specified, using avw.hdr.hist.orient value.\n');
    end
    IMGorient = double(avw.hdr.hist.orient);
end

if ~isfinite(IMGorient),
    if verbose,
        fprintf('...IMGorient is not finite!\n');
    end
    IMGorient = 99;
end

switch IMGorient,
    
case 0, % transverse/axial unflipped
    
    % For the 'transverse unflipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'y' axis                  - from patient posterior to anterior
    % Slices in 'z' axis                  - from patient inferior to superior
    
    if verbose, fprintf('...writing axial unflipped\n'); end
    
    avw.hdr.hist.orient = uint8(0);
    
    SliceDim = double(avw.hdr.dime.dim(4)); % z
    RowDim   = double(avw.hdr.dime.dim(3)); % y
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(4));
    RowSz    = double(avw.hdr.dime.pixdim(3));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    x = 1:PixelDim;
    for z = 1:SliceDim,
        for y = 1:RowDim,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 1, % coronal unflipped
    
    % For the 'coronal unflipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'z' axis                  - from patient inferior to superior
    % Slices in 'y' axis                  - from patient posterior to anterior
    
    if verbose, fprintf('...writing coronal unflipped\n'); end
    
    avw.hdr.hist.orient = uint8(1);
    
    SliceDim = double(avw.hdr.dime.dim(3)); % y
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(3));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    x = 1:PixelDim;
    for y = 1:SliceDim,
        for z = 1:RowDim,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 2, % sagittal unflipped
    
    % For the 'sagittal unflipped' type, the voxels are stored with
    % Pixels in 'y' axis (varies fastest) - from patient posterior to anterior
    % Rows in   'z' axis                  - from patient inferior to superior
    % Slices in 'x' axis                  - from patient right to left
    
    if verbose, fprintf('...writing sagittal unflipped\n'); end
    
    avw.hdr.hist.orient = uint8(2);
    
    SliceDim = double(avw.hdr.dime.dim(2)); % x
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(3)); % y
    SliceSz  = double(avw.hdr.dime.pixdim(2));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(3));
    
    y = 1:PixelDim;
    for x = 1:SliceDim,
        for z = 1:RowDim,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 3, % transverse/axial flipped
    
    % For the 'transverse flipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'y' axis                  - from patient anterior to posterior*
    % Slices in 'z' axis                  - from patient inferior to superior
    
    if verbose,
        fprintf('...writing axial flipped (+Y from Anterior to Posterior)\n');
    end
    
    avw.hdr.hist.orient = uint8(3);
    
    SliceDim = double(avw.hdr.dime.dim(4)); % z
    RowDim   = double(avw.hdr.dime.dim(3)); % y
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(4));
    RowSz    = double(avw.hdr.dime.pixdim(3));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    x = 1:PixelDim;
    for z = 1:SliceDim,
        for y = RowDim:-1:1, % flipped in Y
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 4, % coronal flipped
    
    % For the 'coronal flipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'z' axis                  - from patient inferior to superior
    % Slices in 'y' axis                  - from patient anterior to posterior
    
    if verbose,
        fprintf('...writing coronal flipped (+Z from Superior to Inferior)\n');
    end
    
    avw.hdr.hist.orient = uint8(4);
    
    SliceDim = double(avw.hdr.dime.dim(3)); % y
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(3));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    x = 1:PixelDim;
    for y = 1:SliceDim,
        for z = RowDim:-1:1,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 5, % sagittal flipped
    
    % For the 'sagittal flipped' type, the voxels are stored with
    % Pixels in 'y' axis (varies fastest) - from patient posterior to anterior
    % Rows in   'z' axis                  - from patient superior to inferior
    % Slices in 'x' axis                  - from patient right to left
    
    if verbose,
        fprintf('...writing sagittal flipped (+Z from Superior to Inferior)\n');
    end
    
    avw.hdr.hist.orient = uint8(5);
    
    SliceDim = double(avw.hdr.dime.dim(2)); % x
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(3)); % y
    SliceSz  = double(avw.hdr.dime.pixdim(2));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(3));
    
    y = 1:PixelDim;
    for x = 1:SliceDim,
        for z = RowDim:-1:1, % superior to inferior
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
otherwise, % transverse/axial unflipped
    
    % For the 'transverse unflipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'y' axis                  - from patient posterior to anterior
    % Slices in 'z' axis                  - from patient inferior to superior
    
    if verbose,
        fprintf('...unknown orientation specified, assuming default axial unflipped\n');
    end
    
    avw.hdr.hist.orient = uint8(0);
    
    SliceDim = double(avw.hdr.dime.dim(4)); % z
    RowDim   = double(avw.hdr.dime.dim(3)); % y
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(4));
    RowSz    = double(avw.hdr.dime.pixdim(3));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    x = 1:PixelDim;
    for z = 1:SliceDim,
        for y = 1:RowDim,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
end

fclose(fid);

% Update the header
avw.hdr.dime.dim(2:4) = int16([PixelDim,RowDim,SliceDim]);
avw.hdr.dime.pixdim(2:4) = single([PixelSz,RowSz,SliceSz]);

return







%----------------------------------------------------------------------------

function write_header(fid,avw,verbose)
    
    header_key(fid,avw.hdr.hk);
    image_dimension(fid,avw.hdr.dime);
    data_history(fid,avw.hdr.hist);
    
    % check the file size is 348 bytes
    fbytes = ftell(fid);
    fclose(fid);
    if ~isequal(fbytes,348),
        msg = sprintf('...file size is not 348 bytes!\n');
        warning(msg);
    end
    
return

%----------------------------------------------------------------------------

function header_key(fid,hk)
    
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
    
    fwrite(fid, hk.sizeof_hdr(1),   'int32');  % must be 348!
    
    data_type = sprintf('%-10s',hk.data_type); % ensure it is 10 chars
    fwrite(fid, hk.data_type(1:10), 'uchar');
    
    db_name   = sprintf('%-18s',hk.db_name);   % ensure it is 18 chars
    fwrite(fid, db_name(1:18),      'uchar');
    
    fwrite(fid, hk.extents(1),      'int32');
    fwrite(fid, hk.session_error(1),'int16');
    
    regular   = sprintf('%1s',hk.regular);     % ensure it is 1 char
    fwrite(fid, regular(1),         'uchar');  % might be uint8
    
    %hkey_un0  = sprintf('%1s',hk.hkey_un0);    % ensure it is 1 char
    %fwrite(fid, hkey_un0(1),        'uchar');
    fwrite(fid, hk.hkey_un0(1),     'uint8');
    
    %    >Would you set hkey_un0 as char or uint8?
    %   Really doesn't make any difference.  As far as anyone here can remember,
    %   this was just to pad to an even byte boundary for that structure.  I guess
    %   I'd suggest setting it to a uint8 value of 0 (i.e, truly zero-valued) so
    %   that it doesn't look like anything important!
    %   Denny <hanson.dennis2@mayo.edu>
    
return

%----------------------------------------------------------------------------

function image_dimension(fid,dime)
    
	%struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
	%       char vox_units[4];               /* 16 + 4          */
	%       char cal_units[8];               /* 20 + 8          */
	%       short int unused1;               /* 28 + 2          */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int dim_un0;               /* 34 + 2          */
	%       float pixdim[8];                 /* 36 + 32         */
	%			/*
	%				pixdim[] specifies the voxel dimensions:
	%				pixdim[1] - voxel width
	%				pixdim[2] - voxel height
	%				pixdim[3] - interslice distance
	%					..etc
	%			*/
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
    
	fwrite(fid, dime.dim(1:8),      'int16');
	fwrite(fid, dime.vox_units(1:4),'uchar');
	fwrite(fid, dime.cal_units(1:8),'uchar');
	fwrite(fid, dime.unused1(1),    'int16');
	fwrite(fid, dime.datatype(1),   'int16');
	fwrite(fid, dime.bitpix(1),     'int16');
	fwrite(fid, dime.dim_un0(1),    'int16');
	fwrite(fid, dime.pixdim(1:8),   'float32');
	fwrite(fid, dime.vox_offset(1), 'float32');
    
    % Ensure compatibility with SPM (according to MRIcro)
    if dime.roi_scale == 0, dime.roi_scale = 0.00392157; end
	fwrite(fid, dime.roi_scale(1),  'float32');
    
	fwrite(fid, dime.funused1(1),   'float32');
	fwrite(fid, dime.funused2(1),   'float32');
	fwrite(fid, dime.cal_max(1),    'float32');
	fwrite(fid, dime.cal_min(1),    'float32');
	fwrite(fid, dime.compressed(1), 'int32');
	fwrite(fid, dime.verified(1),   'int32');
	fwrite(fid, dime.glmax(1),      'int32');
	fwrite(fid, dime.glmin(1),      'int32');
	
return

%----------------------------------------------------------------------------

function data_history(fid,hist)
    
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
	
    descrip     = sprintf('%-80s', hist.descrip);       % 80 chars
    aux_file    = sprintf('%-24s', hist.aux_file);      % 24 chars
    originator  = sprintf('%-10s', hist.originator);    % 10 chars
    generated   = sprintf('%-10s', hist.generated);     % 10 chars
    scannum     = sprintf('%-10s', hist.scannum);       % 10 chars
    patient_id  = sprintf('%-10s', hist.patient_id);    % 10 chars
    exp_date    = sprintf('%-10s', hist.exp_date);      % 10 chars
    exp_time    = sprintf('%-10s', hist.exp_time);      % 10 chars
    hist_un0    = sprintf( '%-3s', hist.hist_un0);      %  3 chars
    
    % ---
    % The following should not be necessary, but I actually
    % found one instance where it was, so this totally anal
    % retentive approach became necessary, despite the
    % apparently elegant solution above to ensuring that variables
    % are the right length.
    
    if length(descrip) < 80,
      paddingN = 80-length(descrip);
      padding = char(repmat(double(' '),1,paddingN));
      descrip = [descrip,padding];
    end
    if length(aux_file) < 24,
      paddingN = 24-length(aux_file);
      padding = char(repmat(double(' '),1,paddingN));
      aux_file = [aux_file,padding];
    end
    if length(originator) < 10,
      paddingN = 10-length(originator);
      padding = char(repmat(double(' '),1,paddingN));
      originator = [originator, padding];
    end
    if length(generated) < 10,
      paddingN = 10-length(generated);
      padding = char(repmat(double(' '),1,paddingN));
      generated = [generated, padding];
    end
    if length(scannum) < 10,
      paddingN = 10-length(scannum);
      padding = char(repmat(double(' '),1,paddingN));
      scannum = [scannum, padding];
    end
    if length(patient_id) < 10,
      paddingN = 10-length(patient_id);
      padding = char(repmat(double(' '),1,paddingN));
      patient_id = [patient_id, padding];
    end
    if length(exp_date) < 10,
      paddingN = 10-length(exp_date);
      padding = char(repmat(double(' '),1,paddingN));
      exp_date = [exp_date, padding];
    end
    if length(exp_time) < 10,
      paddingN = 10-length(exp_time);
      padding = char(repmat(double(' '),1,paddingN));
      exp_time = [exp_time, padding];
    end
    if length(hist_un0) < 10,
      paddingN = 10-length(hist_un0);
      padding = char(repmat(double(' '),1,paddingN));
      hist_un0 = [hist_un0, padding];
    end
    
    % -- if you thought that was anal, try this;
    % -- lets check for unusual ASCII char values!
    
    if find(double(descrip)>128),
      indexStrangeChar = find(double(descrip)>128);
      descrip(indexStrangeChar) = ' ';
    end
    if find(double(aux_file)>128),
      indexStrangeChar = find(double(aux_file)>128);
      aux_file(indexStrangeChar) = ' ';
    end
    if find(double(originator)>128),
      indexStrangeChar = find(double(originator)>128);
      originator(indexStrangeChar) = ' ';
    end
    if find(double(generated)>128),
      indexStrangeChar = find(double(generated)>128);
      generated(indexStrangeChar) = ' ';
    end
    if find(double(scannum)>128),
      indexStrangeChar = find(double(scannum)>128);
      scannum(indexStrangeChar) = ' ';
    end
    if find(double(patient_id)>128),
      indexStrangeChar = find(double(patient_id)>128);
      patient_id(indexStrangeChar) = ' ';
    end
    if find(double(exp_date)>128),
      indexStrangeChar = find(double(exp_date)>128);
      exp_date(indexStrangeChar) = ' ';
    end
    if find(double(exp_time)>128),
      indexStrangeChar = find(double(exp_time)>128);
      exp_time(indexStrangeChar) = ' ';
    end
    if find(double(hist_un0)>128),
      indexStrangeChar = find(double(hist_un0)>128);
      hist_un0(indexStrangeChar) = ' ';
    end
    
    
    % --- finally, we write the fields
    
    fwrite(fid, descrip(1:80),    'uchar');
    fwrite(fid, aux_file(1:24),   'uchar');
    
    
    %orient      = sprintf(  '%1s', hist.orient);        %  1 char
    %fwrite(fid, orient(1),        'uchar');
    fwrite(fid, hist.orient(1),   'uint8');     % see note below on char
    
    fwrite(fid, originator(1:10), 'uchar');
    fwrite(fid, generated(1:10),  'uchar');
    fwrite(fid, scannum(1:10),    'uchar');
    fwrite(fid, patient_id(1:10), 'uchar');
    fwrite(fid, exp_date(1:10),   'uchar');
    fwrite(fid, exp_time(1:10),   'uchar');
    fwrite(fid, hist_un0(1:3),    'uchar');
    
    fwrite(fid, hist.views(1),      'int32');
    fwrite(fid, hist.vols_added(1), 'int32');
    fwrite(fid, hist.start_field(1),'int32');
    fwrite(fid, hist.field_skip(1), 'int32');
    fwrite(fid, hist.omax(1),       'int32');
    fwrite(fid, hist.omin(1),       'int32');
    fwrite(fid, hist.smax(1),       'int32');
    fwrite(fid, hist.smin(1),       'int32');
    
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

% See other notes in avw_hdr_read
