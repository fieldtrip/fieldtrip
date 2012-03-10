function [IFF]=openiff(fid,LEN)
% OPENIFF is an auxillary function to SOPEN for 
% opening of IFF files 
% 
% Use SOPEN instead of OPENIFF  
% 
% See also: fopen, SOPEN
%
% References: 

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.

%	$Id$
%	Copyright (C) 2004,2005,2007 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


FLAG.CLOSE = 0; 
if ischar(fid)
        fid = fopen(fid,'r');
        FLAG.CLOSE = 1; 
end;

if nargin<2,
        LEN = 8;
end;		

IFF = [];
K   = 0; 
K1  = 0; 

[tmp,c] = fread(fid,[1,4],'uint8');
while ((LEN>0) | isnan(LEN)) & (c>0),	
        tag     = char(tmp);
        tagsize = fread(fid,1,'uint32');        % which size 
        tagsize0= tagsize + rem(tagsize,2); 
        filepos = ftell(fid);
        LEN = LEN - 8; 
        
%       fprintf(1,'tag: %6s\tpos: %8i\tsize: %8i\n',tag,filepos,tagsize);
        
        if 0,
                VAL = openiff(fid,tagsize);
        elseif strcmp(tag,'FORM')
                [tmp,c] = fread(fid,[1,4],'uint8');tag,
                VAL = setfield([],char(tmp),openiff(fid,tagsize-4));
        elseif 0,strcmp(tag,'RIFF')
                VAL = openiff(fid,tagsize);
        elseif strcmp(tag,'RIFF')
                [tmp,c] = fread(fid,[1,4],'uint8');
                %val = fread(fid,tagsize-4,'uint8');
                val = openiff(fid,tagsize-4);
                VAL = setfield([],char(tmp),val);
        elseif strcmp(tag,'MThd')
                VAL.MThd = fread(fid,tagsize,'uchar');
                VAL.MIDI = openiff(fid,NaN);
                %LEN = NaN;
                %VAL.MIDI = tmp;
 
        elseif strcmp(tag,'LIST');	% CNT_RIFF (EEP 3.1)
                [tmp,c] = fread(fid,[1,4],'uint8');
                VAL = setfield([],char(tmp),openiff(fid,tagsize-4));

        elseif strcmp(tag,'chan')	% CNT_RIFF (EEP 3.1)
                VAL = fread(fid,[1,tagsize/2],'uint16');
        elseif strcmp(tag,'info')	% CNT_RIFF (EEP 3.1)
                VAL = fread(fid,[1,tagsize],'*char');
        elseif strcmp(tag,'eeph')	% CNT_RIFF (EEP 3.1)
                VAL = fread(fid,[1,tagsize],'*char');
        elseif strncmp(tag,'ep  ',4),	% CNT_RIFF (EEP 3.1) 
                VAL = fread(fid,[1,tagsize/4],'uint32');
        elseif strcmp(tag,'hdrl')
                [tmp,c] = fread(fid,[1,4],'uint8');
                VAL = setfield([],char(tmp),openiff(fid,tagsize-4));

        elseif strcmp(tag,'CAT ')
                VAL = openiff(fid,tagsize);

        elseif strcmp(tag,'data')
        %        LEN = fread(fid,1,'uint32');
                VAL = fread(fid,[1,LEN],'*uint8');
                %[tmp,c] = fread(fid,[1,4],'*char');
                %LEN = fread(fid,1,'int32');
                %VAL.data = fread(fid,LEN,'uint8');
                %VAL = openiff(fid,tagsize);
        elseif strncmp(tag,'(c)',3)
                tag = 'Copyright';
                VAL = char(fread(fid,tagsize,'uchar')');
                %VAL.CopyRight = char(VAL);
        else
                if 0,tagsize<1024*2,
                        VAL = fread(fid,tagsize,'uchar');
                else
			VAL = [];
                        VAL.handle = ftell(fid);
                        VAL.size = tagsize; 
                        status = fseek(fid,tagsize,'cof');
                        if status, return; end;
                end
        end;
	        
        if strcmp(tag(3:4),'dc')
                K = K+1;
                ix = (tag(1:2)-48)*[16;1]+1;
                IFF.dc{K,ix} = VAL;
        elseif strcmp(tag(3:4),'wb')
                K1 = K1+1;
                ix = (tag(1:2)-48)*[16;1]+1;
                IFF.wb{K1,ix} = VAL;
        else
                %try,
	                IFF = setfield(IFF,deblank(tag),VAL);
                %catch,
                %end;
        end;
        status = fseek(fid,filepos+tagsize0,'bof');
        LEN = LEN - tagsize0;
        [tmp,c] = fread(fid,[1,4],'uint8');
end;

if ~isfield(IFF,'MThd'), % do not check MIDI files
        if ((LEN ~= 0) & c), 
                fprintf(2,'Warning OPENIFF: LEN=%i %i %i %i %s  %i\n',LEN,filepos,tagsize0,ftell(fid),char(tmp),c);
        end;	
end;

if FLAG.CLOSE,
        fclose(fid);
end;
