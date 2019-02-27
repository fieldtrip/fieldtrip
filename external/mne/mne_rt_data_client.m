classdef mne_rt_data_client < mne_rt_client
    %MNE_RT_DATA_CLIENT is a class to parse incomming fiff data tags send 
    % by mne_rt_server using the fiff data port
    %
    %   This class is inherited of mne_rt_client.
    %
    %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
    %   License : BSD 3-clause
    %
    properties
        m_clientID = -1; % The client ID assigned by mne_rt_server
    end
    
    methods
        % =================================================================
        %% mne_rt_data_client
        function obj = mne_rt_data_client(host, port, numOfRetries)
            %
            % obj = mne_rt_data_client(host, port, numOfRetries)
            %
            % Constructor
            %
            % host          - ip adress of the mne_rt_Server
            % port          - fiff data port of the mne_rt_Server
            % numOfRetries  - number of connection retries, when a
            %                 connection fails
            %
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            if (nargin < 3)
                numOfRetries = 20; % set to -1 for infinite
            end
            obj = obj@mne_rt_client(host, port, numOfRetries); %Superclass call
            obj.getClientId();
        end % mne_rt_data_client
        
        % =================================================================
        %% readInfo
        function [info] = readInfo(obj)
            % reads the measurement info
            import java.net.Socket
            import java.io.*
            
            global FIFF;
            if isempty(FIFF)
                FIFF = fiff_define_constants();
            end
            
            if ~isempty(obj.m_DataInputStream)
                
                t_bReadMeasBlockStart = false;
                t_bReadMeasBlockEnd = false;
                %
                % Find the start
                %
                while(~t_bReadMeasBlockStart)
                    tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);

                    if(tag.kind == FIFF.FIFF_BLOCK_START && tag.data == FIFF.FIFFB_MEAS_INFO)
                        disp('FIFF_BLOCK_START FIFFB_MEAS_INFO'); 
                        t_bReadMeasBlockStart = true;
                    end
                end

                %
                % Parse until the endblock
                %
                
                info.dev_head_t = [];
                info.ctf_head_t = [];
                info.dev_ctf_t = [];
                info.acq_pars = [];
                info.acq_stim = [];
                info.chs = [];
  
                dev_head_t_read = false;
                ctf_head_t_read = false;

                while(~t_bReadMeasBlockEnd)
                    tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                    %
                    %  megacq parameters
                    %
                    if(tag.kind == FIFF.FIFF_BLOCK_START && tag.data == FIFF.FIFFB_DACQ_PARS)                        
                        while(tag.kind ~= FIFF.FIFF_BLOCK_END || tag.data ~= FIFF.FIFFB_DACQ_PARS)
                            tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                            if(tag.kind == FIFF.FIFF_DACQ_PARS)
                                info.acq_pars = tag.data;
                            elseif(tag.kind == FIFF.FIFF_DACQ_STIM)
                                info.acq_stim = tag.data;
                            end
                        end
                    end
                    %
                    %    Coordinate transformations if the HPI result block was not there
                    %
                    if (tag.kind == FIFF.FIFF_COORD_TRANS)
                        if (~dev_head_t_read)
                            info.dev_head_t = tag.data;
                            dev_head_t_read = true;
                        elseif (~ctf_head_t_read)
                            info.ctf_head_t = tag.data;
                            ctf_head_t_read = true;
                        end
                    end
                    %
                    %    Polhemus data
                    %
                    if(tag.kind == FIFF.FIFF_BLOCK_START && tag.data == FIFF.FIFFB_ISOTRAK)
                       info.dig = [];
                       while(tag.kind ~= FIFF.FIFF_BLOCK_END || tag.data ~= FIFF.FIFFB_ISOTRAK)
                            tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                            if(tag.kind == FIFF.FIFF_DIG_POINT)
                                info.dig = [info.dig tag.data];
                            end
                        end
                    end
                    %
                    %    Projectors
                    %
                    if(tag.kind == FIFF.FIFF_BLOCK_START && tag.data == FIFF.FIFFB_PROJ)
                       info.projs = [];
                       while(tag.kind ~= FIFF.FIFF_BLOCK_END || tag.data ~= FIFF.FIFFB_PROJ)
                            tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                            if(tag.kind == FIFF.FIFF_BLOCK_START && tag.data == FIFF.FIFFB_PROJ_ITEM)
                               proj = [];
                               while(tag.kind ~= FIFF.FIFF_BLOCK_END || tag.data ~= FIFF.FIFFB_PROJ_ITEM)
                                    tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                                    switch tag.kind
                                        case FIFF.FIFF_NAME
                                            proj.desc = tag.data;
                                        case FIFF.FIFF_PROJ_ITEM_KIND
                                            proj.kind = tag.data;
                                        case FIFF.FIFF_NCHAN
                                            proj.data.ncol = tag.data;
                                        case FIFF.FIFF_PROJ_ITEM_NVEC
                                            proj.data.nrow = tag.data;
                                        case FIFF.FIFF_MNE_PROJ_ITEM_ACTIVE
                                            proj.active = tag.data;
                                        case FIFF.FIFF_PROJ_ITEM_CH_NAME_LIST
                                            proj.data.col_names = fiff_split_name_list(tag.data);
                                        case FIFF.FIFF_PROJ_ITEM_VECTORS
                                            proj.data.data = tag.data;
                                    end
                                end
                            end
                            
                            if ~isempty(proj)
                                info.projs = [info.projs proj];
                            end
                        end
                    end
                    %
                    %    CTF compensation info
                    %
                    if(tag.kind == FIFF.FIFF_BLOCK_START && tag.data == FIFF.FIFFB_MNE_CTF_COMP)
                       info.comps = [];
                       while(tag.kind ~= FIFF.FIFF_BLOCK_END || tag.data ~= FIFF.FIFFB_MNE_CTF_COMP)
                            tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                            if(tag.kind == FIFF.FIFF_BLOCK_START && tag.data == FIFF.FIFFB_MNE_CTF_COMP_DATA)
                               comp = [];
                               while(tag.kind ~= FIFF.FIFF_BLOCK_END || tag.data ~= FIFF.FIFFB_MNE_CTF_COMP_DATA)
                                    tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                                    switch tag.kind
                                        case FIFF.FIFF_MNE_CTF_COMP_KIND
                                            comp.ctfkind = tag.data;
                                        case FIFF.FIFF_MNE_CTF_COMP_CALIBRATED
                                            comp.save_calibrated = tag.data;
                                        case FIFF.FIFF_MNE_CTF_COMP_DATA
                                            comp.data = tag.data;
                                    end
                                end
                            end
                            
                            if ~isempty(comp)
                                info.comps = [info.comps comp];
                            end
                        end
                    end
                    %
                    %    Bad channels
                    %
                    if(tag.kind == FIFF.FIFF_BLOCK_START && tag.data == FIFF.FIFFB_MNE_BAD_CHANNELS)
                       info.bads = [];
                       while(tag.kind ~= FIFF.FIFF_BLOCK_END || tag.data ~= FIFF.FIFFB_MNE_BAD_CHANNELS)
                            tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                            if(tag.kind == FIFF.FIFF_MNE_CH_NAME_LIST)
                                info.bads = fiff_split_name_list(tag.data);
                            end
                        end
                    end
                    %
                    %    General
                    %
                    if (tag.kind == FIFF.FIFF_SFREQ)
                        info.sfreq = tag.data;
                    elseif (tag.kind == FIFF.FIFF_HIGHPASS)
                        info.highpass = tag.data;
                    elseif (tag.kind == FIFF.FIFF_LOWPASS)
                        info.lowpass = tag.data;
                    elseif (tag.kind == FIFF.FIFF_NCHAN)
                        info.nchan = tag.data;
                    elseif (tag.kind == FIFF.FIFF_MEAS_DATE)
                        info.highpass = tag.data;
                    end
                        
                    
                    if (tag.kind == FIFF.FIFF_CH_INFO)
                        info.chs = [info.chs tag.data];
                    end
                    
                    % END MEAS
                    if(tag.kind == FIFF.FIFF_BLOCK_END && tag.data == FIFF.FIFFB_MEAS_INFO)
                        disp('FIFF_BLOCK_END FIFFB_MEAS_INFO'); 
                        t_bReadMeasBlockEnd = true;
                    end
                end
            else
                error('mne_rt_data_client: no available TcpSocket, call init to establish a connection.');
            end
        end
        
        % =================================================================
        %% readRawBuffer
        function [kind, data] = readRawBuffer(obj, p_nChannels)
            %
            % [kind, data] = readRawBuffer(obj, p_nChannels)
            %
            % reads a raw buffer
            %
            % p_nChannels - number of channels to reshape the incomming 
            %               float array
            %
            % kind        - FIFF_DATA_BUFFER ->
            %               FIFF_BLOCK_START -> data = FIFFB_RAW_DATA
            %               FIFF_BLOCK_END -> data = FIFFB_RAW_DATA
            % data        - the read buffer
            
            %
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %

            import java.net.Socket
            import java.io.*
              
            global FIFF;
            if isempty(FIFF)
                FIFF = fiff_define_constants();
            end
            
            data = [];
            kind = [];
            
            if ~isempty(obj.m_DataInputStream)

                bytes_available = obj.m_DataInputStream.available;
                % Wait for incomming bytes
                while(bytes_available == 0)
                    bytes_available = obj.m_DataInputStream.available;
                end

                tag = mne_rt_data_client.read_tag(obj.m_DataInputStream);
                
                kind = tag.kind;
                
                if(tag.kind == FIFF.FIFF_DATA_BUFFER)
                    nSamples = length(tag.data)/p_nChannels;
                    data = reshape(tag.data, p_nChannels, nSamples);
                else
                    data = tag.data;
                end
            else
                error('mne_rt_data_client: no available TcpSocket, call init to establish a connection.');
            end
        end
        
        % =================================================================
        %% setClientAlias
        function setClientAlias(obj, alias)
            %
            % setClientAlias(obj, alias)
            %
            % sets the alias of the fiff data client -> for convinient
            % identification
            %
            % alias - name which should be used
            %

            %
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %

            import java.net.Socket
            import java.io.*

            global MNE_RT;
            if isempty(MNE_RT)
                MNE_RT = mne_rt_define_commands();
            end
            
            if ~isempty(obj.m_DataOutputStream)
                mne_rt_data_client.sendFiffCommand(obj.m_DataOutputStream, MNE_RT.MNE_RT_SET_CLIENT_ALIAS, alias)
            end
        end
        
        % =================================================================
        %% getClientId
        function [id] = getClientId(obj)
            %
            % [id] = getClientId(obj)
            %
            % sets and returns the client id which was assigned by
            % mne_rt_server

            %
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            
            import java.net.Socket
            import java.io.*
            if(obj.m_clientID == -1)
                global FIFF;
                if isempty(FIFF)
                    FIFF = fiff_define_constants();
                end
                global MNE_RT;
                if isempty(MNE_RT)
                    MNE_RT = mne_rt_define_commands();
                end

                if ~isempty(obj.m_DataOutputStream)

                    mne_rt_data_client.sendFiffCommand(obj.m_DataOutputStream, MNE_RT.MNE_RT_GET_CLIENT_ID)

                    % ID is send as answer
                    tag = obj.read_tag(obj.m_DataInputStream);
                    if (tag.kind == FIFF.FIFF_MNE_RT_CLIENT_ID)
                        obj.m_clientID = tag.data;
                    end                
                end
            end
            id = obj.m_clientID;
        end
    end
    
    methods(Static)
        % =================================================================
        %% sendFiffCommand
        function sendFiffCommand(p_DataOutputStream, p_Cmd, p_data)
            %
            % sendFiffCommand(p_DataOutputStream, p_Cmd, p_data)
            %
            % sends a fiff encoded command to mne_rt_server using the fiff
            % data port
            % command and data are encoded into the fiff tag data
            %
            % p_DataOutputStream - open output data stream
            % p_Cmd              - command to send
            % p_data             - data which should be also encoded
            %

            %
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %

            import java.net.Socket
            import java.io.*
            
            global FIFF;
            if isempty(FIFF)
                FIFF = fiff_define_constants();
            end
            
            if (nargin == 3)
                data = char(p_data);
            elseif(nargin == 2)
                data = [];
            else
                error('Wrong number of arguments.');
            end
            
            kind = FIFF.FIFF_MNE_RT_COMMAND;
            type = FIFF.FIFFT_VOID;
            size = 4+length(data); % first 4 bytes are the command code
            next = 0;
            
            p_DataOutputStream.writeInt(kind);
            p_DataOutputStream.writeInt(type);
            p_DataOutputStream.writeInt(size);
            p_DataOutputStream.writeInt(next);
            p_DataOutputStream.writeInt(p_Cmd); % first 4 bytes are the command code
            if(~isempty(data))
                p_DataOutputStream.writeBytes(data);
            end
            p_DataOutputStream.flush;
        end
        
        % =================================================================
        %% read_tag
        function [tag] = read_tag(p_DataInputStream)
            %
            % [tag] = read_tag(p_DataInputStream)
            %
            % reads a fiff encoded real-time data stream
            %
            % p_DataInputStream - open data stream
            %
            
            %
            %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
            %   License : BSD 3-clause
            %
            import java.net.Socket
            import java.io.*
            
            me='MNE_RT_DATA_CLIENT:read_tag';
            
            %
            % read the tag info
            %
            tagInfo = mne_rt_data_client.read_tag_info(p_DataInputStream);

            %
            % read the tag data
            %
            tag = mne_rt_data_client.read_tag_data(p_DataInputStream, tagInfo);
        end
        
        % =================================================================
        %% read_tag_data
        function [tag] = read_tag_data(p_DataInputStream, p_tagInfo, pos)
        %
        % [tag] = read_tag_data(p_dInputStream, pos)
        %
        % Reads the tag data from a fif stream.
        % if pos is not provided, reading starts from the current stream position
        %
        % p_DataInputStream - the open data stream
        % p_tagInfo         - the tag info
        % pos               - number of bytes to be skipped
        %
        
        %
        %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
        %   License : BSD 3-clause
        %
            import java.net.Socket
            import java.io.*
            
            global FIFF;
            if isempty(FIFF)
                FIFF = fiff_define_constants();
            end

            me='MNE_RT_DATA_CLIENT:read_tag_data';

            if nargin == 3
                p_DataInputStream.skipBytes(pos);
            elseif nargin ~= 2
                error(me,'Incorrect number of arguments');
            end

            tag = p_tagInfo;

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
                        
% Check can't be done in real-time --> moved to the end for reshape
%                         pos = ftell(fid);
%                         fseek(fid,tag.size-4,'cof');
%                         ndim = fread(fid,1,'int32');
%                         fseek(fid,-(ndim+1)*4,'cof');
%                         dims = fread(fid,ndim,'int32');
%                         %
%                         % Back to where the data start
%                         %
%                         fseek(fid,pos,'bof');

                        matrix_type = bitand(data_type,tag.type);

                        el_size = tag.size - 3*4; % 3*4 --> case 2D matrix; ToDo calculate el_size through

                        switch matrix_type
%                             case FIFF.FIFFT_INT
%                                 tag.data = zeros(el_size/4, 1);
%                                 for i = 1:el_size/4
%                                     tag.data(i) = p_DataInputStream.readInt; %idata = fread(fid,dims(1)*dims(2),'int32=>int32');
%                                 end
%                             case FIFF.FIFFT_JULIAN
%                                 tag.data = zeros(el_size/4, 1);
%                                 for i = 1:el_size/4
%                                     tag.data(i) = p_DataInputStream.readInt; %idata = fread(fid,dims(1)*dims(2),'int32=>int32');
%                                 end
                            case FIFF.FIFFT_FLOAT
                                t_MNERTBufferReader = mne_rt_buffer_reader(p_DataInputStream);
                                tmp = t_MNERTBufferReader.readBuffer(el_size);
                                tag.data = typecast(tmp, 'single'); %fdata = fread(fid,dims(1)*dims(2),'single=>double');
                                tag.data = swapbytes(tag.data);
                            otherwise
                                error(me,'Cannot handle a matrix of type %d yet',matrix_type)
                        end
                        
                        % ToDo consider 3D case --> do that by using tag->size
                        dims = zeros(1, 2);
                       
                        dims(1) = p_DataInputStream.readInt;
                        dims(2) = p_DataInputStream.readInt;
                        
                        ndim = p_DataInputStream.readInt;
                        
                        tag.data = reshape(tag.data,dims(1),dims(2))';
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
                        case FIFF.FIFFT_INT
%                             tag.data = zeros(tag.size/4, 1);
%                             for i = 1:tag.size/4
%                                 tag.data(i) = p_DataInputStream.readInt; %fread(fid,tag.size/4,'int32=>int32');
%                             end
                            t_MNERTBufferReader = mne_rt_buffer_reader(p_DataInputStream);
                            tmp = t_MNERTBufferReader.readBuffer(tag.size);
                            tag.data = typecast(tmp, 'int32'); %fread(fid,tag.size/4,'int32=>int32');
                            tag.data = swapbytes(tag.data);
                        case FIFF.FIFFT_FLOAT
                            t_MNERTBufferReader = mne_rt_buffer_reader(p_DataInputStream);
                            tmp = t_MNERTBufferReader.readBuffer(tag.size);
                            tag.data = typecast(tmp, 'single'); %fread(fid,tag.size/4,'single=>double');
                            tag.data = swapbytes(tag.data);
                        case FIFF.FIFFT_STRING
                            t_MNERTBufferReader = mne_rt_buffer_reader(p_DataInputStream);
                            tag.data = t_MNERTBufferReader.readBuffer(tag.size); %fread(fid,tag.size,'uint8=>char')';
                            tag.data = char(tag.data);                            
                        case FIFF.FIFFT_ID_STRUCT
                            tag.data.version = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.machid = zeros(2,1);
                            tag.data.machid(1)  = p_DataInputStream.readInt; %fread(fid,2,'int32=>int32');
                            tag.data.machid(2)  = p_DataInputStream.readInt;
                            tag.data.secs    = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.usecs   = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                        case FIFF.FIFFT_DIG_POINT_STRUCT
                            tag.data.kind    = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.ident   = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.r = zeros(3,1);
                            for i = 1:3
                                tag.data.r(i)	= p_DataInputStream.readFloat; %fread(fid,3,'single=>single');
                            end
                            tag.data.coord_frame = 0;
                        case FIFF.FIFFT_COORD_TRANS_STRUCT
                            tag.data.from = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.to   = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            rot = zeros(9,1);
                            for i = 1:9
                                rot(i) = p_DataInputStream.readFloat; %fread(fid,9,'single=>double');
                            end
                            rot = reshape(rot,3,3)';
                            move = zeros(3,1);
                            for i = 1:3
                                move(i) = p_DataInputStream.readFloat; %fread(fid,3,'single=>double');
                            end
                            tag.data.trans = [ rot move ; [ 0  0 0 1 ]];
                            %
                            % Skip over the inverse transformation
                            % It is easier to just use inverse of trans in Matlab
                            %
                            for i = 1:12 %fseek(fid,12*4,'cof');
                                p_DataInputStream.readFloat;
                            end
                        case FIFF.FIFFT_CH_INFO_STRUCT
                            tag.data.scanno    = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.logno     = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.kind      = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.range     = p_DataInputStream.readFloat; %fread(fid,1,'single=>double');
                            tag.data.cal       = p_DataInputStream.readFloat; %fread(fid,1,'single=>double');
                            tag.data.coil_type = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            %
                            %   Read the coil coordinate system definition
                            %
                            tag.data.loc = zeros(12,1);
                            for i = 1:12
                                tag.data.loc(i) = p_DataInputStream.readFloat; %fread(fid,12,'single=>double');
                            end
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
                            tag.data.unit     = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            tag.data.unit_mul = p_DataInputStream.readInt; %fread(fid,1,'int32=>int32');
                            %
                            %   Handle the channel name
                            %
                            ch_name = zeros(1, 16, 'uint8');
                            for i = 1:16
                                ch_name(i) = p_DataInputStream.readByte;
                            end
                            ch_name   = char(ch_name);
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
                        otherwise
                            error(me,'Unimplemented tag data type %d',tag.type);

                    end
                end
            end

            % if tag.next ~= FIFF.FIFFV_NEXT_SEQ
            %     fseek(fid,tag.next,'bof');
            % end

            return;
        end
            
        % =================================================================
        %% read_tag_info
        function [tag] = read_tag_info(p_DataInputStream, pos)
        %
        % [tag] = read_tag_info(p_inStream, pos)
        %
        % Read tag info from a fif stream.
        % if pos is not provided, reading starts from the current stream position
        %
        % p_DataInputStream - the open data stream
        % pos               - number of bytes to be skipped
        %
        
        %
        %   Author : Christoph Dinh, Matti Hamalainen, MGH Martinos Center
        %   License : BSD 3-clause
        %
            import java.net.Socket
            import java.io.*
            
            global FIFF;
            if isempty(FIFF)
                FIFF = fiff_define_constants();
            end

            me='MNE_RT_DATA_CLIENT:read_tag_info';
            
            if nargin == 2
                p_DataInputStream.skipBytes(pos);
            elseif nargin ~= 1
                error(me,'Incorrect number of arguments');
            end
            
%             while true
%                 bytes_available = p_inStream.available;
%                 if(bytes_available >= 16)
            tag.kind = p_DataInputStream.readInt;
            tag.type = p_DataInputStream.readInt;
            tag.size = p_DataInputStream.readInt;
            tag.next = p_DataInputStream.readInt;
%                     break;
%                 end
%                 % pause 100ms before retrying
%                 pause(0.1);
%             end
            
            return;
        end                    
    end
end

