function [tag] = fiff_read_tag(fid,pos)
%
% [tag] = fiff_read_tag(fid,pos)
%
% Read one tag from a fif file.
% if pos is not provided, reading starts from the current file position
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.16  2008/11/17 15:07:00  msh
%   Added reading of short, unsigned short, and unsigned int data
%
%   Revision 1.15  2008/05/07 10:39:13  msh
%   Added FIFFT_JULIAN
%
%   Revision 1.14  2008/03/25 18:13:29  msh
%   Added reading of sparse matrices
%
%   Revision 1.13  2006/09/24 18:52:43  msh
%   Added FIFFV_REF_MEG_CH to fiff_define_constants.
%   Added coord_frame field to dig point structure.
%
%   Revision 1.12  2006/09/23 13:31:55  msh
%   Added reading of complex and double complex arrays and matrices
%
%   Revision 1.11  2006/05/03 19:03:19  msh
%   Eliminated the use of cast function for Matlab 6.5 compatibility
%
%   Revision 1.10  2006/05/03 18:53:04  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.9  2006/04/26 19:50:58  msh
%   Added fiff_read_mri
%
%   Revision 1.8  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.7  2006/04/17 11:52:15  msh
%   Added coil definition stuff
%
%   Revision 1.6  2006/04/15 12:21:00  msh
%   Several small improvements
%
%   Revision 1.5  2006/04/13 21:20:06  msh
%   Added new routines for the channel information transformation.
%
%   Revision 1.4  2006/04/13 17:53:31  msh
%   Added coordinate frame to the channel info structure
%
%   Revision 1.3  2006/04/13 17:49:12  msh
%   Added some more comments
%
%   Revision 1.2  2006/04/13 17:47:50  msh
%   Make coil coordinate transformation or EEG electrode location out of the channel  info.
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%
global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_tag';

if nargin == 2
    fseek(fid,pos,'bof');
elseif nargin ~= 1
    error(me,'Incorrect number of arguments');
end

tag.kind = fread(fid,1,'int32');
tag.type = fread(fid,1,'uint32');
tag.size = fread(fid,1,'int32');
tag.next = fread(fid,1,'int32');
%
%   The magic hexadecimal values
%
is_matrix           = 4294901760; % ffff0000
matrix_coding_dense = 16384;      % 4000
matrix_coding_CCS   = 16400;      % 4010
matrix_coding_RCS   = 16416;      % 4020
data_type           = 65535;      % ffff
%
if tag.size > 0
    matrix_coding = bitand(is_matrix,tag.type);
    if matrix_coding ~= 0
        matrix_coding = bitshift(matrix_coding,-16);
        %
        %   Matrices
        %
        if matrix_coding == matrix_coding_dense
            %
            % Find dimensions and return to the beginning of tag data
            %
            pos = ftell(fid);
            fseek(fid,tag.size-4,'cof');
            ndim = fread(fid,1,'int32');
            fseek(fid,-(ndim+1)*4,'cof');
            dims = fread(fid,ndim,'int32');
            %
            % Back to where the data start
            %
            fseek(fid,pos,'bof');
            
            matrix_type = bitand(data_type,tag.type);
            
            if ndim == 2
                switch matrix_type
                    case FIFF.FIFFT_INT
                        idata = fread(fid,dims(1)*dims(2),'int32=>int32');
                        tag.data = reshape(idata,dims(1),dims(2))';
                    case FIFF.FIFFT_JULIAN
                        idata = fread(fid,dims(1)*dims(2),'int32=>int32');
                        tag.data = reshape(idata,dims(1),dims(2))';
                    case FIFF.FIFFT_FLOAT
                        fdata = fread(fid,dims(1)*dims(2),'single=>double');
                        tag.data = reshape(fdata,dims(1),dims(2))';
                    case FIFF.FIFFT_DOUBLE
                        ddata = fread(fid,dims(1)*dims(2),'double=>double');
                        tag.data = reshape(ddata,dims(1),dims(2))';
                    case FIFF.FIFFT_COMPLEX_FLOAT
                        fdata = fread(fid,2*dims(1)*dims(2),'single=>double');
                        nel = length(fdata);
                        fdata = complex(fdata(1:2:nel),fdata(2:2:nel));
                        %
                        %   Note: we need the non-conjugate transpose here
                        %
                        tag.data = transpose(reshape(fdata,dims(1),dims(2)));
                    case FIFF.FIFFT_COMPLEX_DOUBLE
                        ddata = fread(fid,2*dims(1)*dims(2),'double=>double');
                        nel = length(ddata);
                        ddata = complex(ddata(1:2:nel),ddata(2:2:nel));
                        %
                        %   Note: we need the non-conjugate transpose here
                        %
                        tag.data = transpose(reshape(ddata,dims(1),dims(2)));
                    otherwise
                        error(me,'Cannot handle a 2D matrix of type %d yet',matrix_type)
                end
            elseif ndim == 3
                switch matrix_type
                    case FIFF.FIFFT_INT
                        idata = fread(fid,dims(1)*dims(2)*dims(3),'int32=>int32');
                        tag.data = reshape(idata,dims(1),dims(2),dims(3));
                    case FIFF.FIFFT_JULIAN
                        idata = fread(fid,dims(1)*dims(2)*dims(3),'int32=>int32');
                        tag.data = reshape(idata,dims(1),dims(2),dims(3));
                    case FIFF.FIFFT_FLOAT
                        fdata = fread(fid,dims(1)*dims(2)*dims(3),'single=>double');
                        tag.data = reshape(fdata,dims(1),dims(2),dims(3));
                    case FIFF.FIFFT_DOUBLE
                        ddata = fread(fid,dims(1)*dims(2)*dims(3),'double=>double');
                        tag.data = reshape(ddata,dims(1),dims(2),dims(3));
                    case FIFF.FIFFT_COMPLEX_FLOAT
                        fdata = fread(fid,2*dims(1)*dims(2)*dims(3),'single=>double');
                        nel = length(fdata);
                        fdata = complex(fdata(1:2:nel),fdata(2:2:nel));
                        tag.data = reshape(fdata,dims(1),dims(2),dims(3));
                    case FIFF.FIFFT_COMPLEX_DOUBLE
                        ddata = fread(fid,2*dims(1)*dims(2)*dims(3),'double=>double');
                        nel = length(ddata);
                        ddata = complex(ddata(1:2:nel),ddata(2:2:nel));
                        tag.data = reshape(ddata,dims(1),dims(2),dims(3));
                    otherwise
                        error(me,'Cannot handle a 3D matrix of type %d yet',matrix_type)
                end
                %
                %   Permute
                %
                tag.data = permute(tag.data,[ 3 2 1 ]);
            else
                error(me, ...
                    'Only two and three dimensional matrices are supported at this time');
            end
        elseif (matrix_coding == matrix_coding_CCS || matrix_coding == matrix_coding_RCS)
            %
            % Find dimensions and return to the beginning of tag data
            %
            pos = ftell(fid);
            fseek(fid,tag.size-4,'cof');
            ndim = fread(fid,1,'int32');
            fseek(fid,-(ndim+2)*4,'cof');
            dims = fread(fid,ndim+1,'int32');
            if ndim ~= 2
                error(me,'Only two-dimensional matrices are supported at this time');
            end
            %
            % Back to where the data start
            %
            fseek(fid,pos,'bof');
            nnz   = dims(1);
            nrow  = dims(2);
            ncol  = dims(3);
            sparse_data = zeros(nnz,3);
            sparse_data(:,3) = fread(fid,nnz,'single=>double');
            if (matrix_coding == matrix_coding_CCS)
                %
                %    CCS
                %
                sparse_data(:,1)  = fread(fid,nnz,'int32=>double') + 1;
                ptrs  = fread(fid,ncol+1,'int32=>double') + 1;
                p = 1;
                for j = 1:ncol
                    while p < ptrs(j+1)
                        sparse_data(p,2) = j;
                        p = p + 1;
                    end
                end
            else
                %
                %    RCS
                %
                sparse_data(:,2)  = fread(fid,nnz,'int32=>double') + 1;
                ptrs  = fread(fid,nrow+1,'int32=>double') + 1;
                p = 1;
                for j = 1:nrow
                    while p < ptrs(j+1)
                        sparse_data(p,1) = j;
                        p = p + 1;
                    end
                end
            end
            tag.data = spconvert(sparse_data);
            tag.data(nrow,ncol) = 0.0;
        else
            error(me,'Cannot handle other than dense or sparse matrices yet')
        end
    else
        %
        %   All other data types
        %
        switch tag.type
            %
            %   Simple types
            %
            case FIFF.FIFFT_BYTE
                tag.data = fread(fid,tag.size,'uint8=>uint8');
            case FIFF.FIFFT_SHORT
                tag.data = fread(fid,tag.size/2,'int16=>int16');
            case FIFF.FIFFT_INT
                tag.data = fread(fid,tag.size/4,'int32=>int32');
            case FIFF.FIFFT_USHORT
                tag.data = fread(fid,tag.size/2,'uint16=>uint16');
            case FIFF.FIFFT_UINT
                tag.data = fread(fid,tag.size/4,'uint32=>uint32');
            case FIFF.FIFFT_FLOAT
                tag.data = fread(fid,tag.size/4,'single=>double');
            case FIFF.FIFFT_DOUBLE
                tag.data = fread(fid,tag.size/8,'double');
            case FIFF.FIFFT_STRING
                tag.data = fread(fid,tag.size,'uint8=>char')';
            case FIFF.FIFFT_DAU_PACK16
                tag.data = fread(fid,tag.size/2,'int16=>int16');
            case FIFF.FIFFT_COMPLEX_FLOAT
                tag.data = fread(fid,tag.size/4,'single=>double');
                nel = length(tag.data);
                tag.data = complex(tag.data(1:2:nel),tag.data(2:2:nel));
            case FIFF.FIFFT_COMPLEX_DOUBLE
                tag.data = fread(fid,tag.size/8,'double');
                nel = length(tag.data);
                tag.data = complex(tag.data(1:2:nel),tag.data(2:2:nel));
                %
                %   Structures
                %
            case FIFF.FIFFT_ID_STRUCT
                tag.data.version = fread(fid,1,'int32=>int32');
                tag.data.machid  = fread(fid,2,'int32=>int32');
                tag.data.secs    = fread(fid,1,'int32=>int32');
                tag.data.usecs   = fread(fid,1,'int32=>int32');
            case FIFF.FIFFT_DIG_POINT_STRUCT
                tag.data.kind    = fread(fid,1,'int32=>int32');
                tag.data.ident   = fread(fid,1,'int32=>int32');
                tag.data.r       = fread(fid,3,'single=>single');
                tag.data.coord_frame = 0;
            case FIFF.FIFFT_COORD_TRANS_STRUCT
                tag.data.from = fread(fid,1,'int32=>int32');
                tag.data.to   = fread(fid,1,'int32=>int32');
                rot  = fread(fid,9,'single=>double');
                rot = reshape(rot,3,3)';
                move = fread(fid,3,'single=>double');
                tag.data.trans = [ rot move ; [ 0  0 0 1 ]];
                %
                % Skip over the inverse transformation
                % It is easier to just use inverse of trans in Matlab
                %
                fseek(fid,12*4,'cof');
            case FIFF.FIFFT_CH_INFO_STRUCT
                tag.data.scanno    = fread(fid,1,'int32=>int32');
                tag.data.logno     = fread(fid,1,'int32=>int32');
                tag.data.kind      = fread(fid,1,'int32=>int32');
                tag.data.range     = fread(fid,1,'single=>double');
                tag.data.cal       = fread(fid,1,'single=>double');
                tag.data.coil_type = fread(fid,1,'int32=>int32');
                %
                %   Read the coil coordinate system definition
                %
                tag.data.loc        = fread(fid,12,'single=>double');
                tag.data.coil_trans  = [];
                tag.data.eeg_loc     = [];
                tag.data.coord_frame = FIFF.FIFFV_COORD_UNKNOWN;
                %
                %   Convert loc into a more useful format
                %
                loc = tag.data.loc;
                if tag.data.kind == FIFF.FIFFV_MEG_CH || tag.data.kind == FIFF.FIFFV_REF_MEG_CH
                    tag.data.coil_trans  = [ [ loc(4:6) loc(7:9) loc(10:12) loc(1:3) ] ; [ 0 0 0 1 ] ];
                    tag.data.coord_frame = FIFF.FIFFV_COORD_DEVICE;
                elseif tag.data.kind == FIFF.FIFFV_EEG_CH
                    if norm(loc(4:6)) > 0
                        tag.data.eeg_loc     = [ loc(1:3) loc(4:6) ];
                    else
                        tag.data.eeg_loc = [ loc(1:3) ];
                    end
                    tag.data.coord_frame = FIFF.FIFFV_COORD_HEAD;
                end
                %
                %   Unit and exponent
                %
                tag.data.unit     = fread(fid,1,'int32=>int32');
                tag.data.unit_mul = fread(fid,1,'int32=>int32');
                %
                %   Handle the channel name
                %
                ch_name   = fread(fid,16,'uint8=>char')';
                %
                % Omit nulls
                %
                len = 16;
                for k = 1:16
                    if ch_name(k) == 0
                        len = k-1;
                        break
                    end
                end
                tag.data.ch_name = ch_name(1:len);
            case FIFF.FIFFT_OLD_PACK
                offset   = fread(fid,1,'single=>double');
                scale    = fread(fid,1,'single=>double');
                tag.data = fread(fid,(tag.size-8)/2,'int16=>short');
                tag.data = scale*single(tag.data) + offset;
            case FIFF.FIFFT_DIR_ENTRY_STRUCT
                tag.data = struct('kind',{},'type',{},'size',{},'pos',{});
                for k = 1:tag.size/16-1
                    kind = fread(fid,1,'int32');
                    type = fread(fid,1,'uint32');
                    tagsize = fread(fid,1,'int32');
                    pos  = fread(fid,1,'int32');
                    tag.data(k).kind = kind;
                    tag.data(k).type = type;
                    tag.data(k).size = tagsize;
                    tag.data(k).pos  = pos;
                end
                
            otherwise
                error(me,'Unimplemented tag data type %d',tag.type);
                
        end
    end
end

if tag.next ~= FIFF.FIFFV_NEXT_SEQ
    fseek(fid,tag.next,'bof');
end

return;

end
