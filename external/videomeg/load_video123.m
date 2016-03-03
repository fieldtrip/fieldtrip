function [ts, ids, site_id, is_sender, frame_sizes, frame_ptrs] = load_video123(fname, outfldr)
% OUTFLDR is the output folder for storing jpeg images of individual
% frames. If OUTFLDR is an emptry string, do not create any jpeg images
    
%--------------------------------------------------------------------------
%   Copyright (C) 2015 BioMag Laboratory, Helsinki University Central Hospital
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, version 3.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
    
    MAGIC_STR = 'ELEKTA_VIDEO_FILE';

    inpf = fopen(fname, 'rb');
    
    % get the file size
    fseek(inpf, 0, 'eof');
    inp_sz = ftell(inpf);
    fseek(inpf, 0, 'bof');
    
    header_str = fread(inpf, length(MAGIC_STR), 'uchar');
    
    % make sure the file is a valid video file
    assert(strcmp(char(header_str'), MAGIC_STR));
    
    ver = fread(inpf, 1, 'uint32');
    assert(ver==1 || ver==2 || ver==3);
    
    if ver==3
        site_id = fread(inpf, 1, 'uint8');
        sender_flag = fread(inpf, 1, 'uint8');
        is_sender = (sender_flag==1);
    else
        site_id = -1;
        is_sender = false;
    end 
    
    ts = [];
    ids = [];
    frame_sizes = [];
    frame_ptrs = [];
    
    while(ftell(inpf) < inp_sz)
        tstamp = fread(inpf, 1, 'uint64=>uint64');
        if(ver > 1)
            chunk_id = fread(inpf, 1, 'uint64=>uint64');
        else
            chunk_id = -1;
        end
                
        chunk_sz = fread(inpf, 1, 'uint32');
        frame_ptrs = [frame_ptrs, ftell(inpf)];
        chunk = fread(inpf, chunk_sz, 'uchar=>uchar');
        
        ts = [ts, tstamp];
        ids = [ids, chunk_id];
        frame_sizes = [frame_sizes, chunk_sz];
        
        if(~isempty(outfldr))
            outf = fopen(sprintf('%s/%06i.jpg', outfldr, length(ts)), 'w+');
            fwrite(outf, chunk, 'uchar');
            fclose(outf);
        end
    end
    
    fclose(inpf);
end
