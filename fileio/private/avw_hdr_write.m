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

% $Revision: 883 $ $Date: 2009/01/14 09:24:45 $

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
    version = '[$Revision: 883 $]';
    fprintf('AVW_HDR_WRITE [v%s]\n',version(12:16));  tic;
end

%----------------------------------------------------------------------------
% Check inputs

if ~exist('avw','var'),
    warning('...no input avw - calling avw_hdr_make\n');
    avw = avw_hdr_make;
elseif isempty(avw),
    warning('...empty input avw - calling avw_hdr_make\n');
    avw = avw_hdr_make;
elseif ~isfield(avw,'hdr'),
    warning('...empty input avw.hdr - calling avw_hdr_make\n');
    avw = avw_hdr_make;
end
if ~isequal(avw.hdr.hk.sizeof_hdr,348),
    msg = sprintf('...avw.hdr.hk.sizeof_hdr must be 348!\n');
    error(msg);
end

quit = 0;
if ~exist('fileprefix','var'),
    if isfield(avw,'fileprefix'),
        if ~isempty(avw.fileprefix),
            fileprefix = avw.fileprefix;
        else,
            quit = 1;
        end
    else
        quit = 1;
    end
    if quit,
        helpwin avw_hdr_write;
        error('...no input fileprefix - see help avw_hdr_write\n\n');
        return;
    end
end

if findstr('.hdr',fileprefix),
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
if fid < 0,
    msg = sprintf('Cannot write to file %s.hdr\n',fileprefix);
    error(msg);
else
    if verbose, fprintf('...writing %s Analyze header.\n',machine); end
    write_header(fid,avw,verbose);
end

if verbose, t=toc; fprintf('...done (%5.2f sec).\n\n',t); end

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
    %           /*
    %               pixdim[] specifies the voxel dimensions:
    %               pixdim[1] - voxel width
    %               pixdim[2] - voxel height
    %               pixdim[3] - interslice distance
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
