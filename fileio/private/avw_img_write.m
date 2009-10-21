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

% $Revision: 1.1 $ $Date: 2009/01/14 09:24:45 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%                    The Analyze format is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------
% Check inputs

if ~exist('avw','var'),
    doc avw_img_write;
    error('...no input avw.');
elseif isempty(avw),
    error('...empty input avw.');
elseif ~isfield(avw,'img'),
    error('...empty input avw.img');
end

if ~exist('fileprefix','var'),
    if isfield(avw,'fileprefix'),
        if ~isempty(avw.fileprefix),
            fileprefix = avw.fileprefix;
        else
            fileprefix = [];
        end
    else
        fileprefix = [];
    end
end
if isempty(fileprefix),
    [fileprefix, pathname, filterindex] = uiputfile('*.hdr','Specify an output Analyze .hdr file');
    if pathname, cd(pathname); end
    if ~fileprefix,
        doc avw_img_write;
        error('no output .hdr file specified');
    end
end

if findstr('.hdr',fileprefix),
%    fprintf('AVW_IMG_WRITE: Removing .hdr extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.hdr','');
end
if findstr('.img',fileprefix),
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
if verbose,
    version = '[$Revision: 1.1 $]';
    fprintf('\nAVW_IMG_WRITE [v%s]\n',version(12:16));  tic;
end

fid = fopen(sprintf('%s.img',fileprefix),'w',machine);
if fid < 0,
    msg = sprintf('Cannot open file %s.img\n',fileprefix);
    error(msg);
else
    avw = write_image(fid,avw,fileprefix,IMGorient,machine,verbose);
end

if verbose,
    t=toc; fprintf('...done (%5.2f sec).\n\n',t);
end

% MUST write header after the image, to ensure any
% orientation changes during image write are saved
% in the header
avw_hdr_write(avw,fileprefix,machine,verbose);

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

switch double(avw.hdr.dime.datatype),
case   1,
    avw.hdr.dime.bitpix = int16( 1); precision = 'bit1';
case   2,
    avw.hdr.dime.bitpix = int16( 8); precision = 'uchar';
case   4,
    avw.hdr.dime.bitpix = int16(16); precision = 'int16';
case   8,
    avw.hdr.dime.bitpix = int16(32); precision = 'int32';
case  16,
    avw.hdr.dime.bitpix = int16(32); precision = 'single';
case  32,
    error('...complex datatype not yet supported.\n');
case  64,
    avw.hdr.dime.bitpix = int16(64); precision = 'double';
case 128,
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
