function [S,HDR,time] = sread(HDR,NoS,StartPos)
% SREAD loads selected segments of signal file
%
% [S,HDR] = sread(HDR [,NoS [,StartPos]] )
% NoS       Number of seconds, default = 1 (second)
% StartPos  Starting position, if not provided the following data is read continously from the file. 
%                    no reposition of file pointer is performed
%
% HDR=sopen(Filename,'r',CHAN);
% [S,HDR] = sread(HDR, NoS, StartPos)
%      	reads NoS seconds beginning at StartPos
% 
% [S,HDR] = sread(HDR, inf) 
%      	reads til the end starting at the current position 
% 
% [S,HDR] = sread(HDR, N*HDR.Dur) 
%	reads N trials of an BKR file 
% 
%
% See also: fread, SREAD, SWRITE, SCLOSE, SSEEK, SREWIND, STELL, SEOF

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.

%	$Id$
%	(C) 1997-2005,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

S = [];
time = []; 

if nargin<2, 
        NoS = inf; 
end;

if ~isnumeric(NoS) || (NoS<0),
        fprintf(HDR.FILE.stderr,'Error SREAD: NoS must be non-negative number\n');
        return;
end;
if (nargin>=3) 
        if (StartPos<0),
                fprintf(HDR.FILE.stderr,'Error SREAD: StartPos must be non-negative\n');
                return;
        end;
        tmp = HDR.SampleRate*StartPos;
        if tmp ~= round(tmp),
        %        fprintf(HDR.FILE.stderr,'Warning SREAD: StartPos yields non-integer position\n');
                StartPos = round(tmp)/HDR.SampleRate;
        end;
else
        StartPos = HDR.FILE.POS/HDR.SampleRate; 
end;

tmp = HDR.SampleRate*NoS;
if tmp ~= round(tmp),
        fprintf(HDR.FILE.stderr,'Warning SREAD: NoS yields non-integer position [%f, %f]\n',NoS,HDR.SampleRate);
        NoS = round(tmp)/HDR.SampleRate;
end;

% define HDR.out.EVENT. This is used by EEGLAB. 
ix = (HDR.EVENT.POS >= StartPos*HDR.SampleRate) & (HDR.EVENT.POS <= (StartPos+NoS)*HDR.SampleRate); 
HDR.out.EVENT.POS = HDR.EVENT.POS(ix)-StartPos;
HDR.out.EVENT.TYP = HDR.EVENT.TYP(ix);
if isfield(HDR.EVENT,'CHN')
        if ~isempty(HDR.EVENT.CHN)
                HDR.out.EVENT.CHN = HDR.EVENT.CHN(ix);
        end;
end;
if isfield(HDR.EVENT,'DUR')
        if ~isempty(HDR.EVENT.DUR)
                HDR.out.EVENT.DUR = HDR.EVENT.DUR(ix);
        end;
end;

STATUS = 0; 
if isfield(HDR,'THRESHOLD')     % save THRESHOLD status (will be modified in BKR);
        THRESHOLD = HDR.THRESHOLD; 
end
FLAG_CALIB_DONE = 0; 

if 0, 
elseif strcmp(HDR.TYPE,'BDF'),
        if nargin==3,
                HDR.FILE.POS = round(HDR.SampleRate*StartPos);
        end;

        nr     = min(HDR.NRec*HDR.SPR-HDR.FILE.POS, NoS*HDR.SampleRate);
	block1 = floor(HDR.FILE.POS/HDR.SPR);
	ix1    = HDR.FILE.POS - block1*HDR.SPR;	% starting sample (minus one) within 1st block 
	nb     = ceil((HDR.FILE.POS+nr)/HDR.SPR)-block1;
    	fp     = HDR.HeadLen + block1*HDR.AS.bpb;
    	status = fseek(HDR.FILE.FID, fp, 'bof');

	count  = 0;
        if HDR.NS==0,

        else
                if (HDR.AS.spb*nb<=2^22), % faster access
			S  = repmat(NaN,HDR.SPR*nb,length(HDR.InChanSelect)); 
                        [s,c] = fread(HDR.FILE.FID,[3*HDR.AS.spb, nb],'uint8');
                        s = reshape(2.^[0,8,16]*reshape(s(:),3,c/3),[HDR.AS.spb, nb]);
                        c = c/3;
                        for k = 1:length(HDR.InChanSelect),
                               K = HDR.InChanSelect(k);
                               if (HDR.AS.SPR(K)>0)
                                       S(:,k) = rs(reshape(s(HDR.AS.bi(K)+1:HDR.AS.bi(K+1),:),HDR.AS.SPR(K)*nb,1),HDR.AS.SPR(K),HDR.SPR);
                               else
                                       S(:,k) = NaN;
                               end;
                        end;
                        S = S(ix1+1:ix1+nr,:);
                        count = nr;

                        if HDR.FLAG.OVERFLOWDETECTION,  % BDF overflow detection is based on Status bit20
	                        K = HDR.BDF.Status.Channel;
        	                OVERFLOW = ~bitand(reshape(s(HDR.AS.bi(K)+1:HDR.AS.bi(K+1),:),HDR.AS.SPR(K)*nb,1),2^19);
        	                OVERFLOW = rs(OVERFLOW,HDR.AS.SPR(K),HDR.SPR);
	                        OVERFLOW = OVERFLOW(ix1+1:ix1+nr,:);
        	                S(OVERFLOW>0,:)=NaN; 
        	        end;        
                else
			S      = repmat(NaN,nr,length(HDR.InChanSelect)); 
                        while (count<nr);
                                len   = ceil(min([(nr-count)/HDR.SPR,2^22/HDR.AS.spb]));
                                [s,c] = fread(HDR.FILE.FID,[3*HDR.AS.spb, len],'uint8=>uint8');
                                s1    = zeros(HDR.SPR*c/(3*HDR.AS.spb),length(HDR.InChanSelect));
                                for k = 1:length(HDR.InChanSelect), 
                                        K = HDR.InChanSelect(k);
                                        tmp = 2.^[0,8,16]*double(reshape(s(HDR.AS.bi(K)*3+1:HDR.AS.bi(K+1)*3,:),3,HDR.AS.SPR(K)*c/HDR.AS.bpb));
                                        if (HDR.AS.SPR(K)>0)
                                                s1(:,k) = rs(tmp',HDR.AS.SPR(K),HDR.SPR);
                                        else
                                                s1(:,k) = NaN;
                                        end;
                                end;
	                        if HDR.FLAG.OVERFLOWDETECTION,  % BDF overflow detection is based on Status bit20
		                        K = HDR.BDF.Status.Channel;
                                        tmp = 2.^[0,8,16]*double(reshape(s(HDR.AS.bi(K)*3+1:HDR.AS.bi(K+1)*3,:),3,HDR.AS.SPR(K)*c/HDR.AS.bpb));
                                        OVERFLOW = rs(~bitand(tmp',2^19),HDR.AS.SPR(K),HDR.SPR);
	        	                s1(OVERFLOW>0,:)=NaN; 
	        	        end;        
                                ix2   = min(nr-count, size(s1,1)-ix1);
                                S(count+1:count+ix2,:) = s1(ix1+1:ix1+ix2,:);
                                count = count+ix2;
                                ix1   = 0; 
                        end;	
                end;
                HDR.FILE.POS = HDR.FILE.POS + count;
                S = S - 2^24*(S>=2^23);
        end

elseif strcmp(HDR.TYPE,'EDF') || strcmp(HDR.TYPE,'GDF') || strcmp(HDR.TYPE,'BDF') || strcmp(HDR.TYPE,'ACQ'),
	% experimental, might replace SDFREAD.M 
        if nargin==3,
                HDR.FILE.POS = round(HDR.SampleRate*StartPos);
        end;

        nr     = min(HDR.NRec*HDR.SPR-HDR.FILE.POS, NoS*HDR.SampleRate);
	S      = repmat(NaN,nr,length(HDR.InChanSelect)); 

	block1 = floor(HDR.FILE.POS/HDR.SPR);
	ix1    = HDR.FILE.POS- block1*HDR.SPR;	% starting sample (minus one) within 1st block 
	nb     = ceil((HDR.FILE.POS+nr)/HDR.SPR)-block1;
    	fp     = HDR.HeadLen + block1*HDR.AS.bpb;
    	STATUS = fseek(HDR.FILE.FID, fp, 'bof');
        count  = 0;
        if HDR.NS==0,
        elseif all(HDR.GDFTYP==HDR.GDFTYP(1)),
                if (HDR.AS.spb*nb<=2^24), % faster access
                        S = [];
                        %[HDR.AS.spb, nb,block1,ix1,HDR.FILE.POS],
                        [s,c] = fread(HDR.FILE.FID,[HDR.AS.spb, nb],gdfdatatype(HDR.GDFTYP(1)));
                        for k = 1:length(HDR.InChanSelect),
                               K = HDR.InChanSelect(k);
                               if (HDR.AS.SPR(K)>0)
                                       S(:,k) = rs(reshape(s(HDR.AS.bi(K)+1:HDR.AS.bi(K+1),:),HDR.AS.SPR(K)*nb,1),HDR.AS.SPR(K),HDR.SPR);
                               else
                                       S(:,k) = NaN;
                               end;
                        end;
                        S = S(ix1+1:ix1+nr,:);
                        count = nr;
                else
                        S = repmat(NaN,[nr,length(HDR.InChanSelect)]);
                        while (count<nr);
                                len   = ceil(min([(nr-count)/HDR.SPR,2^22/HDR.AS.spb]));
                                [s,c] = fread(HDR.FILE.FID,[HDR.AS.spb, len],gdfdatatype(HDR.GDFTYP(1)));
                                s1    = zeros(HDR.SPR*c/HDR.AS.spb,length(HDR.InChanSelect));
                                for k = 1:length(HDR.InChanSelect), 
                                        K = HDR.InChanSelect(k);
                                        if (HDR.AS.SPR(K)>0)
                                                tmp = reshape(s(HDR.AS.bi(K)+1:HDR.AS.bi(K+1),:),HDR.AS.SPR(K)*c/HDR.AS.spb,1);
                                                s1(:,k) = rs(tmp,HDR.AS.SPR(K),HDR.SPR);
                                        end;
                                end;
                                ix2   = min(nr-count, size(s1,1)-ix1);
                                S(count+1:count+ix2,:) = s1(ix1+1:ix1+ix2,:);
                                count = count+ix2;
                                ix1   = 0; 
                        end;	
                end;
        else
                fprintf(2,'SREAD (GDF): different datatypes - this might take some time.\n');
                
                S = repmat(NaN,[nr,length(HDR.InChanSelect)]);
                while (count<nr);
                        s = [];
                        for k=1:length(HDR.AS.TYP),
                                [s0,tmp] = fread(HDR.FILE.FID,[HDR.AS.c(k), 1],gdfdatatype(HDR.AS.TYP(k)));
                                s = [s;s0];
                        end;
                        
                        s1    = repmat(NaN,[HDR.SPR,length(HDR.InChanSelect)]);
                        for k = 1:length(HDR.InChanSelect), 
                                K = HDR.InChanSelect(k);
                                if (HDR.AS.SPR(K)>0)
                                        s1(:,k) = rs(s(HDR.AS.bi(K)+1:HDR.AS.bi(K+1),:),HDR.AS.SPR(K),HDR.SPR);
                                end;
                        end;
                        ix2   = min(nr-count, size(s1,1)-ix1);
                        S(count+1:count+ix2,:) = s1(ix1+1:ix1+ix2,:);
                        count = count+HDR.SPR;
                        ix1   = 0; 
                end;	
        end;
        if strcmp(HDR.TYPE,'GDF')       % read non-equidistant sampling channels of GDF2.0 format
                if (HDR.VERSION>1.94) %& isfield(HDR.EVENT,'VAL'),
                        for k = 1:length(HDR.InChanSelect), 
                                ch = HDR.InChanSelect(k);
                                if (HDR.AS.SPR(ch)==0),
                                        ix = find((HDR.EVENT.TYP==hex2dec('7fff')) & (HDR.EVENT.CHN==ch));
                                        pix= HDR.EVENT.POS(ix)-HDR.FILE.POS;
                                        ix1= find((pix > 0) & (pix <= count));
                                        S(pix(ix1),k)=HDR.EVENT.DUR(ix(ix1));
                                end;
                        end;
                end;
        end
        HDR.FILE.POS = HDR.FILE.POS + count;
        
        
elseif strcmp(HDR.TYPE,'AINF'),
        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.SampleRate*StartPos*2*(HDR.NS+2),'bof');        
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;

        %[S,count] = fread(HDR.FILE.FID,[HDR.NS+2,HDR.SampleRate*NoS],'int16');
        nr = min(HDR.SampleRate*NoS, HDR.SPR*HDR.NRec-HDR.FILE.POS);
        S  = []; 
	time = [];
        count = 0; 
        while (count<nr),
               	[s,c] = fread(HDR.FILE.FID, [HDR.NS+2, min(nr-count,floor(2^24/HDR.NS))], 'int16');
		if nargout>2,
			time  = [time; [s(1:2,:)'+2^16*(s(1:2,:)'<0)]*(2.^[16;0])];
		end;	
               	S = [S; s(2+HDR.InChanSelect,:)'];
              	count = count + c/(HDR.NS+2); 
	end; 
        HDR.FILE.POS = HDR.FILE.POS + count;
        
        
elseif strmatch(HDR.TYPE,{'BKR'}),
        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.SampleRate*HDR.NS*StartPos*2,'bof');        
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        [S,count] = fread(HDR.FILE.FID,[HDR.NS,HDR.SampleRate*NoS],'int16');
        if ~isempty(S),
                S = S(HDR.InChanSelect,:)';
                HDR.FILE.POS = HDR.FILE.POS + count/HDR.NS;
	else
		S = zeros(0,length(HDR.InChanSelect));		
        end;
        THRESHOLD(HDR.AS.TRIGCHAN,:)=NaN; % do not apply overflow detection for Trigger channel 

        
elseif strmatch(HDR.TYPE,{'AIF','SND','WAV','Sigma'})
        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.SampleRate*HDR.AS.bpb*StartPos,'bof');
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;

        maxsamples = min(HDR.SPR*HDR.NRec - HDR.FILE.POS, HDR.SampleRate*NoS);
        [S,count] = fread(HDR.FILE.FID,[HDR.NS,maxsamples], gdfdatatype(HDR.GDFTYP));

        S = S(HDR.InChanSelect,:)';
        HDR.FILE.POS = HDR.FILE.POS + count/HDR.NS;
        
        if ~HDR.FLAG.UCAL,
                if isfield(HDR.FILE,'TYPE')
                        if HDR.FILE.TYPE==1,
                                S = mu2lin(S);
                        end;
                end;
        end;

        
elseif strmatch(HDR.TYPE,{'BLSC2','CFWB','CNT','DEMG','DDT','ET-MEG','ISHNE','Nicolet','RG64'}),

	tc = strcmp(HDR.TYPE,'CFWB') && isfield(HDR,'FLAG') && isfield(HDR.FLAG,'TimeChannel') && HDR.FLAG.TimeChannel;

        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.AS.bpb*round(HDR.SampleRate*StartPos),'bof');
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        if strcmpi(HDR.FLAG.OUTPUT,'single'),
		DT = [gdfdatatype(HDR.GDFTYP),'=>single'];
        elseif any(HDR.GDFTYP==[1:6,16]) && ~exist('OCTAVE_VERSION','builtin'),
        	% preserve data type
		DT = ['*',gdfdatatype(HDR.GDFTYP)];
	else
		% convert to double
		DT = [gdfdatatype(HDR.GDFTYP)];
        end;
        maxsamples = min(HDR.SampleRate*NoS, HDR.NRec*HDR.SPR-HDR.FILE.POS);
	S = []; count = 0;
	while maxsamples>0,
		% the maximum block size of 2^23 is a heuristical value 
    		[s,c] = fread(HDR.FILE.FID, [HDR.NS+tc,min(2^23/HDR.NS,maxsamples)], DT);
    		c = c/(HDR.NS+tc);
		count = count + c;
		maxsamples = maxsamples - c;
        	if c>0,
            		S = [S; s(HDR.InChanSelect+tc,:)'];
            	else 
		        fprintf(HDR.FILE.stderr,'Warning SREAD(%s): could not read %i samples, only %i samples read\n',HDR.TYPE,maxsamples+count,count);
            		break; 	
    		end;
        end;
	HDR.FILE.POS = HDR.FILE.POS + count;


elseif strcmp(HDR.TYPE,'EPL'),
        if nargin==3,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        startblock = floor(HDR.FILE.POS/HDR.SPR);
        STATUS   = fseek(HDR.FILE.FID, HDR.HeadLen + startblock*HDR.AS.bpb, 'bof'); % fseek needed because HDR.FILE.POS can be changed by SSEEK
        curblock = startblock;
        endpos   = min(HDR.FILE.POS+NoS*HDR.SampleRate, HDR.NRec*HDR.SPR);

        [datablock,count] = fread(HDR.FILE.FID, [(1+HDR.NS)*HDR.SPR,ceil(endpos/HDR.SPR)-startblock], 'int16');
        datablock = reshape(datablock(256+1:end,:),HDR.NS,HDR.SPR*size(datablock,2))'; % remove mark track, and reshape data  

        S = datablock(HDR.FILE.POS-startblock*HDR.SPR+1:endpos-startblock*HDR.SPR,HDR.InChanSelect); 
        HDR.FILE.POS = HDR.FILE.POS + size(S,1);


elseif strcmp(HDR.TYPE,'SMA'),
        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.SampleRate*HDR.AS.bpb*StartPos,'bof');        
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        tmp = min(NoS*HDR.SampleRate,(HDR.NRec*HDR.SPR-HDR.FILE.POS));
        [S,count] = fread(HDR.FILE.FID,[HDR.NS,tmp],'float32'); % read data frame
        tmp = HDR.NS*tmp;
        if count < tmp,
                fprintf(HDR.FILE.stderr,'Warning SREAD SMA: only %i out of %i samples read\n',count/HDR.NS,tmp/HDR.NS);
        end;
        S = S(HDR.InChanSelect,:)';
        HDR.FILE.POS = HDR.FILE.POS + count/HDR.NS;
        
        HDR.SMA.events = diff(sign([HDR.Filter.T0',S(HDR.SMA.EVENT_CHANNEL,:)]-HDR.SMA.EVENT_THRESH))>0;
        HDR.EVENT.POS = find(HDR.SMA.events);
        HDR.EVENT.TYP = HDR.SMA.events(HDR.EVENT.POS);
        
        if size(S,2) > 0,
                HDR.Filter.T0 = S(HDR.SMA.EVENT_CHANNEL,size(S,2))';
        end;
        
        
elseif strcmp(HDR.TYPE,'RDF'),
        S = [];
        if nargin>2,
                HDR.FILE.POS = StartPos;
        end;
        POS = HDR.FILE.POS;
        
        NoSeg = min(NoS,length(HDR.Block.Pos)-HDR.FILE.POS);
        count = 0;
        S = zeros(NoSeg*HDR.SPR, length(HDR.InChanSelect));
        
        for k = 1:NoSeg,
                STATUS = fseek(HDR.FILE.FID,HDR.Block.Pos(POS+k),-1);
                
                % Read nchans and block length
                tmp = fread(HDR.FILE.FID,34+220,'uint16');
                
                %STATUS = fseek(HDR.FILE.FID,2,0);
                nchans = tmp(2); %fread(HDR.FILE.FID,1,'uint16');
                %fread(HDR.FILE.FID,1,'uint16');
                block_size = tmp(4); %fread(HDR.FILE.FID,1,'uint16');
                %ndupsamp = fread(HDR.FILE.FID,1,'uint16');
                %nrun = fread(HDR.FILE.FID,1,'uint16');
                %err_detect = fread(HDR.FILE.FID,1,'uint16');
                %nlost = fread(HDR.FILE.FID,1,'uint16');
                nevents = tmp(9); %fread(HDR.FILE.FID,1,'uint16');
                %STATUS = fseek(HDR.FILE.FID,50,0);
                
                [data,c] = fread(HDR.FILE.FID,[nchans,block_size],'int16');
                %S = [S; data(HDR.InChanSelect,:)']; 	% concatenate data blocks
                S((k-1)*HDR.SPR+(1:c/nchans),:) = data(HDR.InChanSelect,:)';
                count = count + c;
        end;
        HDR.FILE.POS = HDR.FILE.POS + NoSeg; 
        
        
elseif strcmp(HDR.TYPE,'LABVIEW'),
        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.SampleRate*HDR.AS.bpb*StartPos,'bof');        
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        [S,count] = fread(HDR.FILE.FID,[HDR.NS,HDR.SampleRate*NoS],'int32');
        if count,
                S = S(HDR.InChanSelect,:)';
                HDR.FILE.POS = HDR.FILE.POS + count/HDR.NS;
        end;
        
        
elseif strcmp(HDR.TYPE,'alpha'),
        if nargin==3,
                POS = HDR.SampleRate*StartPos*HDR.AS.bpb/HDR.SPR;
                if POS~=ceil(POS),
                        fprintf(HDR.FILE.stderr,'warning SREAD (alpha): starting position is non-integer (%f)\n',POS);     
                end
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen + POS,'bof');        
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        
        nr = min(HDR.SampleRate*NoS, HDR.NRec*HDR.SPR-HDR.FILE.POS)*HDR.AS.bpb;
        if (nr - round(nr)) < .01,
    		nr = round(nr);
	else
	        fprintf(HDR.FILE.stderr,'Error SREAD (alpha): can not deal with odd number of samples \n');     
                return;
        end
        
        if HDR.Bits==12,
                [s,count] = fread(HDR.FILE.FID,[3,nr/3],'uint8');
                s(1,:) = s(1,:)*16 + floor(s(2,:)/16); 	
                s(3,:) = s(3,:)+ mod(s(2,:),16)*256; 	
                s = reshape(s([1,3],:),2*size(s,2),1);
                s = s - (s>=2^11)*2^12;
		nr = floor(length(s)/HDR.NS);
                S = reshape(s(1:nr*HDR.NS),HDR.NS,nr);
                count = count*2/3;
                
        elseif HDR.Bits==16,
                [S,count] = fread(HDR.FILE.FID,[HDR.NS,nr],'int16');
                
        elseif HDR.Bits==32,
                [S,count] = fread(HDR.FILE.FID,[HDR.NS,nr],'int32');
        end;        
        
        if count,
                S = S(HDR.InChanSelect,:)';
                HDR.FILE.POS = HDR.FILE.POS + count/HDR.NS;
        end;

                
elseif strcmp(HDR.TYPE,'MIT'),
        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.SampleRate*HDR.AS.bpb*StartPos,'bof');        
                tmp = HDR.SampleRate*StartPos;
                if HDR.FILE.POS~=tmp,
                        HDR.mode8.accu = zeros(1,length(HDR.InChanSelect));
                        HDR.mode8.valid= 0;
                end;
        end;

        DataLen = NoS*HDR.SampleRate/HDR.SPR;
        if HDR.VERSION == 212, 
                [A,count] = fread(HDR.FILE.FID, [1,DataLen*HDR.AS.bpb], 'uint8');  % matrix with 3 rows, each 8 bits long, = 2*12bit
		DataLen = floor(count/HDR.AS.bpb);
		if (count~= DataLen*HDR.AS.bpb) && isfinite(DataLen),
			fprintf(HDR.FILE.stderr,'Warning SREAD (MIT): non-integer block length %i,%f\n',count, DataLen*HDR.AS.bpb);
			%HDR = sseek(HDR,HDR.FILE.POS,'bof');
			%return;
                        A = A(1:DataLen*HDR.AS.bpb);
		end;
                S = [A(1:3:end) + mod(A(2:3:end),16)*256; A(3:3:end) + floor(A(2:3:end)/16)*256]; 
                clear A;
                S = S - 2^12*(S>=2^11);	% 2-th complement
                S = reshape(S,HDR.AS.spb,prod(size(S))/HDR.AS.spb)';
                
        elseif HDR.VERSION == 310, 
                [A,count] = fread(HDR.FILE.FID, [HDR.AS.bpb/2, DataLen], 'uint16'); 
                A = A'; DataLen = count/HDR.AS.bpb*2; 
                for k = 1:ceil(HDR.AS.spb/3),
                        k1=3*k-2; k2=3*k-1; k3=3*k;
                        S(:,3*k-2) = floor(mod(A(:,k*2-1),2^12)/2);	
                        S(:,3*k-1) = floor(mod(A(:,k*2),2^12)/2);	
                        S(:,3*k  ) = floor(A(:,k*2-1)*(2^-11)) + floor(A(:,k*2)*(2^-11))*2^5; 
                        S = mod(S(:,1:HDR.AS.spb),2^10);
                        S = S - 2^10*(S>=2^9);	% 2-th complement
                end;
		S = S;
                
        elseif HDR.VERSION == 311, 
                [A,count] = fread(HDR.FILE.FID, [HDR.AS.bpb/4, DataLen], 'uint32');
                A = A'; DataLen = count/HDR.AS.bpb*4;
                for k = 1:ceil(HDR.AS.spb/3),
                        S(:,3*k-2) = mod(A(:,k),2^10);	
                        S(:,3*k-1) = mod(floor(A(:,k)*2^(-11)),2^10);	
                        S(:,3*k)   = mod(floor(A(:,k)*2^(-22)),2^10);	
                        S = S(:,1:HDR.AS.spb);
                        S = S - 2^10*(S>=2^9);	% 2-th complement
                end;
		S = S';
                
        elseif HDR.VERSION == 8, 
                [S,count] = fread(HDR.FILE.FID, [HDR.AS.spb,DataLen], 'int8');  
                S = S'; DataLen = count/HDR.AS.spb          
                
        elseif HDR.VERSION == 80, 
                [S,count] = fread(HDR.FILE.FID, [HDR.AS.spb,DataLen], 'uint8');  
                S = S'-128; DataLen = count/HDR.AS.spb;
                
        elseif HDR.VERSION == 160, 
                [S,count] = fread(HDR.FILE.FID, [HDR.AS.spb,DataLen], 'uint16');  
                S = S'-2^15; DataLen = count/HDR.AS.spb;
                
        elseif HDR.VERSION == 16, 
                [S,count] = fread(HDR.FILE.FID, [HDR.AS.spb,DataLen], 'int16'); 
                S = S'; DataLen = count/HDR.AS.spb;
                
        elseif HDR.VERSION == 61, 
                [S,count] = fread(HDR.FILE.FID, [HDR.AS.spb,DataLen], 'int16'); 
                S = S'; DataLen = count/HDR.AS.spb;
                
        else
                fprintf(2, 'ERROR MIT-ECG: format %i not supported.\n',HDR.VERSION); 
                
        end;
        if any(HDR.AS.SPR>1),
                A = S;
                S = zeros(size(A,1)*HDR.SPR,length(HDR.InChanSelect));
                for k = 1:length(HDR.InChanSelect),
                        ch = HDR.InChanSelect(k);
                        ix = HDR.AS.bi(ch)+1:HDR.AS.bi(ch+1);
                        S(:,k)=rs(reshape(A(:,ix)',size(A,1)*HDR.AS.SPR(ch),1),HDR.AS.SPR(ch),HDR.SPR);
                end
        else
                S = S(:,HDR.InChanSelect);
        end
        if HDR.VERSION == 8, 
                if HDR.FILE.POS==0,
                        HDR.mode8.accu = zeros(1,length(HDR.InChanSelect));
                        HDR.mode8.valid= 1;
                end; 
                if ~HDR.mode8.valid;
                        fprintf(2,'Warning SREAD: unknown offset (TYPE=MIT, mode=8) \n');
                else
                        S(1,:) = S(1,:) + HDR.mode8.accu;
                end;        
                S = cumsum(S);
                HDR.mode8.accu = S(size(S,1),:);
	end;
        HDR.FILE.POS = HDR.FILE.POS + DataLen;   	
        
        
elseif strcmp(HDR.TYPE,'TMS32'),
        if nargin==3,
                HDR.FILE.POS = round(HDR.SampleRate*StartPos);
        end;

	blockN = floor(HDR.FILE.POS/HDR.SPR);
	ix1    = HDR.FILE.POS- blockN*HDR.SPR;	% starting sample (minus one) within 1st block 
    	fp     = HDR.HeadLen + blockN*HDR.AS.bpb;
    	status = fseek(HDR.FILE.FID, fp, 'bof');

        nr     = min(HDR.AS.endpos-HDR.FILE.POS, NoS*HDR.SampleRate);
	S      = repmat(NaN,nr,length(HDR.InChanSelect)); 
	count  = 0;

        fread(HDR.FILE.FID,86,'uint8');
    	while (count<nr) && ~feof(HDR.FILE.FID),
                if all(HDR.GDFTYP==HDR.GDFTYP(1))
                        [s,c] = fread(HDR.FILE.FID,[HDR.NS,HDR.SPR],gdfdatatype(HDR.GDFTYP(1)));
                else
                        s = repmat(NaN,HDR.NS,HDR.SPR);
			c = 0;
                        for k1 = 1:HDR.SPR,
                                for k2 = 1:HDR.NS,
                                        [s(k2,k1),c2] = fread(HDR.FILE.FID,1,gdfdatatype(HDR.GDFTYP(k2)));
                                end;
                        end;
                end;
		ix2 = min(nr-count, size(s,2)-ix1);
		S(count+1:count+ix2,:) = s(HDR.InChanSelect, ix1+1:ix1+ix2)';
		count = count + ix2; 
		ix1 = 0;	% reset starting index, 
                fread(HDR.FILE.FID,86,'uint8');
        end;
	HDR.FILE.POS = HDR.FILE.POS + count;
        
        
elseif strcmp(HDR.TYPE,'EGI'),
        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.AS.bpb*StartPos,'bof');        
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        
        if HDR.FLAG.TRIGGERED,
                NoS = min(NoS,(HDR.NRec-HDR.FILE.POS));
                S = zeros(NoS*HDR.SPR,length(HDR.InChanSelect))+NaN;
                for i = (1:NoS),
                        SegmentCatIndex(HDR.FILE.POS+i) = fread(HDR.FILE.FID,1,'uint16');
                        SegmentStartTime(HDR.FILE.POS+i) = fread(HDR.FILE.FID,1,'uint32');
                        
                        [s,count] = fread(HDR.FILE.FID, [HDR.NS + HDR.EGI.N, HDR.SPR*HDR.NRec], gdfdatatype(HDR.GDFTYP));
                        tmp = (HDR.NS + HDR.EGI.N) * HDR.SPR;
	                if isfinite(tmp) && (count < tmp),
                                fprintf(HDR.FILE.stderr,'Warning SREAD EGI: only %i out of %i samples read\n',count,tmp);
                        end;
                        HDR.FILE.POS = HDR.FILE.POS + count/tmp;
                        
                        if (HDR.EGI.N > 0),
                                [HDR.EVENT.POS,HDR.EVENT.CHN,HDR.EVENT.TYP] = find(s(HDR.NS+1:size(s,1),:)');
	                        HDR.EVENT.DUR = ones(size(HDR.EVENT.POS)); 
                        end 
                        S((i-1)*HDR.SPR + (1:size(s,2)),:) = s(HDR.InChanSelect,:)';
                end;
        else
                [S,count] = fread(HDR.FILE.FID,[HDR.NS + HDR.EGI.N, HDR.SampleRate*NoS],gdfdatatype(HDR.GDFTYP));
                tmp = HDR.SampleRate * NoS;
                if isfinite(tmp) && (count < tmp),
                        fprintf(HDR.FILE.stderr,'Warning SREAD EGI: only %i out of %i samples read\n',count,tmp);
                end;
                HDR.FILE.POS = HDR.FILE.POS + round(count/(HDR.NS + HDR.EGI.N));
                S = S(HDR.InChanSelect,:)';
        end;
        
        
elseif strcmp(HDR.TYPE,'AVG'),
        S = repmat(nan,HDR.SPR,HDR.NS);
        count = 0;
        for i = 1:HDR.NS, 
                [tmp,c]     = fread(HDR.FILE.FID,5,'uint8'); % no longer used 
                count = count + c;
                [S(:,i), c] = fread(HDR.FILE.FID,HDR.SPR,'float');
                count = count + c*4;
        end
        S = S(:,HDR.InChanSelect);
        HDR.FILE.POS = HDR.FILE.POS + count/HDR.AS.bpb;
        
        
elseif strcmp(HDR.TYPE,'COH'),
        warning('.COH data not tested yet')
        if (prod(size(NoS))==1) && (nargin>2), 
                rows = NoS; cols = StartPos;
        elseif prod(size(NoS))==2
                rows = NoS(1); cols = NoS(2);
        else
                fprintf(HDR.FILE.stderr,'Error SREAD mode=COH: invalid arguments.\n');
        end;
        
        STATUS = fseek(HDR.FILE.FID,HDR.COH.directory(rows,cols)+8,'bof'); % skip over a small unused header of 8 bytes 
        sr = fread(HDR.FILE.FID, HDR.SPR, 'float32');  % read real part of coherence    
        si = fread(HDR.FILE.FID, HDR.SPR, 'float32');  % read imag part of coherence    
        S = sr + i * si;
        
        
elseif strcmp(HDR.TYPE,'CSA'),
        warning('.CSA data not tested yet')
        S = fread(HDR.FILE.FID, [HDR.NRec*(HDR.SPR+6)*HDR.NS], 'float32');	        
        
        
elseif strcmp(HDR.TYPE,'EEG'),
        if nargin>2,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.AS.bpb*StartPos,'bof');        
        end;
        
        NoS = min(NoS, HDR.NRec-HDR.FILE.POS);
        S   = zeros(NoS*HDR.SPR, length(HDR.InChanSelect));
        count = 0;
        for i = 1:NoS, %h.compsweeps,
                h.sweep(i).accept   = fread(HDR.FILE.FID,1,'uchar');
                tmp		    = fread(HDR.FILE.FID,2,'ushort');
                h.sweep(i).ttype    = tmp(1);
                h.sweep(i).correct  = tmp(2);
                h.sweep(i).rt       = fread(HDR.FILE.FID,1,'float32');
                tmp  		    = fread(HDR.FILE.FID,2,'ushort');
                h.sweep(i).response = tmp(1);
                h.sweep(i).reserved = tmp(2);
                
                [signal,c] = fread(HDR.FILE.FID, [HDR.NS,HDR.SPR], gdfdatatype(HDR.GDFTYP));
                
                S(i*HDR.SPR+(1-HDR.SPR:0),:) = signal(HDR.InChanSelect,:)';
                count = count + c;
        end;
        HDR.FILE.POS = HDR.FILE.POS + count/HDR.AS.spb;        
        
        
elseif strcmp(HDR.TYPE,'MFER'),
	if (HDR.FRAME.N ~= 1),
		fprintf(2,'Warning MWFOPEN: files with more than one frame not implemented, yet.\n');
		return;
	end
	
	N = 1;
	if ~isfield(HDR,'data'),
		STATUS = fseek(HDR.FILE.FID,HDR.FRAME.POS(N),'bof');
		[tmp,count] = fread(HDR.FILE.FID,HDR.FRAME.sz(N,1:2),gdfdatatype(HDR.FRAME.TYP(N)));
        	if isnan(HDR.NRec),
        		HDR.NRec = count/(HDR.SPR*HDR.NS);
        	end;

        	if count==(HDR.SPR*HDR.NS), %% alternate mode format
        		tmp = reshape(tmp,[HDR.SPR,HDR.NS]);
        	else
        	        tmp = reshape(tmp,[HDR.SPR,HDR.NS,HDR.NRec]);   % convert into 3-Dim
        	        tmp = permute(tmp,[1,3,2]);                     % re-order dimensions
        	        tmp = reshape(tmp,[HDR.SPR*HDR.NRec,HDR.NS]);   % make 2-Dim 
        	end;
		HDR.data = tmp;
	end;

	if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        
        nr = min(HDR.SampleRate*NoS,size(HDR.data,1)-HDR.FILE.POS);
	S  = HDR.data(HDR.FILE.POS + (1:nr), HDR.InChanSelect);
        HDR.FILE.POS = HDR.FILE.POS + nr;
	
        
elseif strcmp(HDR.TYPE,'BCI2000'),
        if nargin==3,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.SampleRate*HDR.AS.bpb*StartPos,'bof');        
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        [S,count] = fread(HDR.FILE.FID,[HDR.NS,HDR.SampleRate*NoS],HDR.BCI2000.GDFTYP,HDR.BCI2000.StateVectorLength);
        if count,
                S = S(HDR.InChanSelect,:)';
                HDR.FILE.POS = HDR.FILE.POS + count/HDR.NS;
        end;

        
elseif strcmp(HDR.TYPE,'native') || strcmp(HDR.TYPE,'SCP'),
	if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;

        nr = min(round(HDR.SampleRate * NoS), size(HDR.data,1) - HDR.FILE.POS);
	S = HDR.data(HDR.FILE.POS+1:HDR.FILE.POS+nr,:);
        HDR.FILE.POS = HDR.FILE.POS + nr;
%	FLAG_CALIB_DONE = 1; 
        
elseif strcmp(HDR.TYPE,'NEX'),
        %% hack: read NEX data once, and transform into "native"
        for k = 1:HDR.NEX.NS,
                fseek(HDR.FILE.FID, HDR.NEX.offset(k), 'bof');
                if 0,  % elseif HDR.NEX.type(k)==1, 
                        
                elseif HDR.NEX.type(k)==2,  % interval
                        tmp = fread(HDR.FILE.FID, [HDR.NEX.nf(k),2], 'int32');
                        HDR.EVENT.POS = [HDR.EVENT.POS; tmp(:,1)];
                        HDR.EVENT.DUR = [HDR.EVENT.DUR; tmp(:,2)];
                        HDR.EVENT.CHN = [HDR.EVENT.CHN; repmat(k,size(tmp,1),1)];
                        
                elseif HDR.NEX.type(k)==3,  % waveform 
                        HDR.HeadLen = HDR.NEX.offset(k);
                        HDR.NEX6.ts{k} = fread(HDR.FILE.FID, [HDR.NEX.nf(k),1], 'int32');
                        HDR.NEX6.data{k} = fread(HDR.FILE.FID, [HDR.NEX.SPR(k), HDR.NEX.nf(k)], 'int16');
                        
                elseif HDR.NEX.type(k)==5, % continous variable  
                        HDR.NEX5.ts{k} = fread(HDR.FILE.FID, [HDR.NEX.nf(k), 2], 'int32')';
                        HDR.data(:,HDR.AS.chanreduce(k)) = fread(HDR.FILE.FID, [HDR.NEX.SPR(k), 1], 'int16');
                        
                elseif HDR.NEX.type(k)==6,  % marker
                        ts = fread(HDR.FILE.FID, [1,HDR.NEX.nf(k)], 'int32');
                        names = zeros(1,64);
                        m = zeros(HDR.NEX.SPR(k), nl, nm);
                        for j=1:nm
                                names(j, :) = fread(HDR.FILE.FID, [1 64], 'uint8');
                                for p = 1:HDR.NEX.SPR(k)
                                        m(p, :, j) = fread(HDR.FILE.FID, [1 nl], 'uint8');
                                end
                        end
                        HDR.NEX.names = names;
                        HDR.NEX.m = m;
                else
                        ts = fread(HDR.FILE.FID, [1,HDR.NEX.nf(k)], 'int32');
                        HDR.NEX0.ts{k}=ts;
                end;
        end
        fclose(HDR.FILE.FID);
        HDR.FILE.OPEN = 0; 
        HDR.TYPE = 'native';
	HDR.data = HDR.data(:,HDR.InChanSelect);
        
        % sequence for reading "native" format
	if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
		
        nr = min([round(HDR.SampleRate * NoS), size(HDR.data,1) - HDR.FILE.POS]);
        S  = HDR.data(HDR.FILE.POS + (1:nr), :);
        HDR.FILE.POS = HDR.FILE.POS + nr;
        
        
elseif strcmp(HDR.TYPE,'PLEXON'),
        %% hack: read PLEXON data once, and transform into "native"
        HDR.data = repmat(NaN,max(HDR.PLX.adcount),HDR.NS);
        NRec = 0; 
        nET  = 0; 
        ET   = repmat(NaN,1024,6);
        wav  = zeros(1024,HDR.PLX.wavlen); 
        
        tscount=zeros(5,130);
        wfcount=zeros(5,130);
        evcount=zeros(1,300);
        adcount=zeros(1,212);
        
        adpos = zeros(1,HDR.NS);
        typ_ubyte = fread(HDR.FILE.FID,2,'int16');
        while ~feof(HDR.FILE.FID),
                TYP  = typ_ubyte(1);
                ubyte= typ_ubyte(2);
                POS  = fread(HDR.FILE.FID,1,'int32');
                tmp  = fread(HDR.FILE.FID,4,'int16');
                CHN  = tmp(1)+1; 
                unit = tmp(2);
                nwf  = tmp(3); 
                DUR  = tmp(4);
    	        if nwf>0,
		        wf = fread(HDR.FILE.FID,[1,DUR],'int16');
                end;
		if TYP<5,
			nET = nET + 1; 
    		        if size(ET,1) < nET;       % memory allocation
        	                ET  = [ET ; repmat(NaN,size(ET))];
        	        end;
        	        ET(nET,:) = [TYP,ubyte*2^32+POS,CHN,DUR,unit,nwf];
		end;
                if TYP==1, % spike
                        tscount(unit+1,CHN) = tscount(unit+1,CHN)+1; 
                        if DUR>0,
                                wfcount(unit+1,CHN) = wfcount(unit+1,CHN)+1; 
                        end;
                elseif TYP==4, % events, 
                        evcount(CHN) = evcount(CHN) + 1; 
                elseif TYP==5,  % continous
                        t1 = adpos(CHN) + 1;
                        t2 = adpos(CHN) + DUR;
                        HDR.data(t1:t2,CHN) = wf';
                        adpos(CHN) = t2;
                end;
                NRec = NRec+1;
                typ_ubyte = fread(HDR.FILE.FID,2,'int16');
        end
        fclose(HDR.FILE.FID);
        HDR.FILE.OPEN = 0; 
        
        HDR.PLX2.tscount = tscount; 
        HDR.PLX2.wfcount = wfcount; 
        HDR.PLX2.evcount = evcount; 
        HDR.PLX2.adcount = adcount; 
        
        ET = ET(1:nET,:);
        HDR.EVENT.ET  = ET; 
        HDR.EVENT.TYP = ET(:,1); 
        HDR.EVENT.POS = ET(:,2); 
        HDR.EVENT.CHN = ET(:,3); 
        HDR.EVENT.DUR = ET(:,4); 
        HDR.EVENT.unit= ET(:,5); 

        HDR.TYPE = 'native';
	HDR.data = HDR.data(:,HDR.InChanSelect);
        
        % sequence for reading "native" format
	if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;

        nr = min(round(HDR.SampleRate * NoS), size(HDR.data,1) - HDR.FILE.POS);
        S  = HDR.data(HDR.FILE.POS + (1:nr), :);
        HDR.FILE.POS = HDR.FILE.POS + nr;
        
        
elseif strcmp(HDR.TYPE,'SCP'),
	% this branch is not in use yet. 
	% its experiemental to improve performance

	% decompress data in first call 
	HT = HDR.Huffman.HT; 
        
	dd = [0:255]';
        ACC = zeros(size(dd));
        c = 0;
        for k2 = 1:8,
                ACC = ACC + (dd>127).*(2^c);
                dd  = mod(dd*2, 256);
                c   = c + 1;
        end;
        
	S2 = []; 
	for k3 = 5:6,
		if (k3==5) && isfield(HDR,'SCP5');
			SCP = HDR.SCP5;
		elseif (k3==6) && isfield(HDR,'SCP6');
			SCP = HDR.SCP6;
		else 
			SCP = []; 
		end;

		if ~isempty(SCP),
                        if ~isfield(HDR,'SCP2'),
				S2 = SCP.data(:,HDR.InChanSelect);                      
        
                        elseif HDR.SCP2.NHT==19999,
                                HuffTab = HDR.Huffman.DHT;
				S2 = zeros(0,length(HDR.InChanSelect));
                                for k0 = 1:length(HDR.InChanSelect), k = HDR.InChanSelect(k0); %HDR.NS,
                                        s2 = SCP.data{k};
                                        s2 = [s2; repmat(0,ceil(max(HDR.SCP2.HT(:,4))/8),1)];
					k1 = 0;	
					l2 = 0; 
					accu = 0;
					c  = 0; 
					x  = [];
					HT = HDR.SCP2.HT(find(HDR.SCP2.HT(:,1)==1),3:7);
					while (l2 < HDR.LeadPos(k,2)),
						while ((c < max(HT(:,2))) && (k1<length(s2)-1));
							k1 = k1 + 1;
							dd = s2(k1);
							accu = accu + ACC(dd+1)*(2^c);
							c = c + 8;
						end;

                                                ixx = 1;
						acc = accu - 2^32*floor(accu*(2^(-32)));   % bitand returns NaN if accu >= 2^32
						while (bitand(acc,2^HT(ixx,1)-1) ~= HT(ixx,5)),
							ixx = ixx + 1;
						end;
                                                
                                                dd = HT(ixx,2) - HT(ixx,1);
						if HT(ixx,3)==0,
							HT = HDR.SCP2.HT(find(HDR.SCP2.HT(:,1)==HT(ixx,5)),3:7);
							fprintf(HDR.FILE.stderr,'Warning SCPOPEN: Switching Huffman Tables is not tested yet.\n');
						elseif (dd==0),
							l2 = l2 + 1;
							x(l2) = HT(ixx,4);
						else %if (HT(ixx,3)>0),
							l2 = l2 + 1;
							
                                                        tmp = floor(accu*(2^(-HT(ixx,1))));       % 
							%tmp = bitshift(accu,-HT(ixx,1));
                                                        tmp = tmp - (2^dd)*floor(tmp*(2^(-dd)));  % 
							%tmp = bitand(tmp,2^dd)
                                                        
                                                        % reverse bit-pattern
                                                        if dd==8,
                                                                tmp = ACC(tmp+1);
                                                        else
                                                                tmp = dec2bin(tmp);
                                                                tmp = [char(repmat('0',1,dd-length(tmp))),tmp];
                                                                tmp = bin2dec(tmp(length(tmp):-1:1));
                                                        end
                                                        x(l2) = tmp-(tmp>=(2^(dd-1)))*(2^dd);
						end;
						accu = floor(accu*2^(-HT(ixx,2)));
						c = c - HT(ixx,2); 
					end;

                                        if k0==1,
                                                S2 = x';
                                        elseif size(x,2)==size(S2,1),
                                                S2(:,k0) = x';
					else
	                                        fprintf(HDR.FILE.stderr,'Error SCPOPEN: Huffman decoding failed (%i) \n',size(x,1));
						return;
                                        end;
				end;
                                
                                
                        elseif (HDR.SCP2.NHT==19999), % alternative decoding algorithm. 
                                HuffTab = HDR.Huffman.DHT;

                                for k0 = 1:length(HDR.InChanSelect), k = HDR.InChanSelect(k0); %HDR.NS,
                                        tmp  = SCP.data{k};
                                        accu = [tmp(4)+256*tmp(3)+65536*tmp(2)+2^24*tmp(1)];
                                        %accu = bitshift(accu,HDR.SCP2.prefix,32);
                                        c  = 0; %HDR.SCP2.prefix;
                                        l  = 4;
                                        l2 = 0;
                                        clear x;
                                        Ntmp = length(tmp);
                                        tmp = [tmp; zeros(4,1)];
                                        while c <= 32, %1:HDR.SPR(k),
                                                ixx = 1;
                                                while (bitand(accu,HDR.Huffman.mask(ixx)) ~= HDR.Huffman.PREFIX(ixx)), 
                                                        ixx = ixx + 1;
                                                end;

                                                if ixx < 18,
                                                        c = c + HDR.Huffman.prefix(ixx);
                                                        %accu  = bitshift(accu, HDR.Huffman.prefix(ixx),32);
                                                        accu  = mod(accu.*(2^HDR.Huffman.prefix(ixx)),2^32);
                                                        l2    = l2 + 1;
                                                        x(l2) = HuffTab(ixx,1);
                                                        
                                                elseif ixx == 18,
                                                        c = c + HDR.Huffman.prefix(ixx) + 8;
                                                        %accu = bitshift(accu, HDR.Huffman.prefix(ixx),32);
                                                        accu  = mod(accu.*(2^HDR.Huffman.prefix(ixx)),2^32);
                                                        l2    = l2 + 1;
                                                        
                                                        acc1  = mod(floor(accu*2^(-24)),256);
                                                        %accu = bitshift(accu, 8, 32);
                                                        accu  = mod(accu*256, 2^32);
                                                        
                                                        x(l2) = acc1-(acc1>=2^7)*2^8;
                                                        acc2  = 0;
                                                        for kk = 1:8,
                                                                acc2 = acc2*2 + mod(acc1,2);
                                                                acc1 = floor(acc1/2);
                                                        end;
                                                        
                                                elseif ixx == 19,
                                                        c = c + HDR.Huffman.prefix(ixx);
                                                        %accu = bitshift(accu, HDR.Huffman.prefix(ixx),32);
                                                        accu  = mod(accu.*(2^HDR.Huffman.prefix(ixx)),2^32);
                                                        l2    = l2 + 1;
                                                        while (c > 7) && (l < Ntmp),
                                                                l = l+1;
                                                                c = c-8;
                                                                accu = accu + tmp(l)*2^c;
                                                        end;
                                                        
                                                        acc1 = mod(floor(accu*2^(-16)),2^16);
                                                        %accu = bitshift(accu, 16, 32);
                                                        accu = mod(accu.*(2^16), 2^32);
                                                        
                                                        x(l2) = acc1-(acc1>=2^15)*2^16;
                                                        acc2 = 0;
                                                        for kk= 1:16,
                                                                acc2 = acc2*2+mod(acc1,2);
                                                                acc1 = floor(acc1/2);
                                                        end;
                                                        %x(l2) = acc2;
                                                        c = c + 16;
                                                end;
                                                
                                                while (c > 7) && (l < Ntmp),
                                                        l = l+1;
                                                        c = c-8;
                                                        accu = accu + tmp(l)*(2^c);
                                                end;
                                        end;

                                        x = x(1:end-1)';
                                        if k==1,
                                                S2=x;
                                        elseif size(x,1)==size(S2,1),
                                                S2(:,k0) = x;
					else
	                                        fprintf(HDR.FILE.stderr,'Error SCPOPEN: Huffman decoding failed (%i) \n',size(x,1));
						return;
                                        end;
                                end;
                        elseif (HDR.SCP2.NHT==1) && (HDR.SCP2.NCT==1) && (HDR.SCP2.prefix==0), 
				S2 = SCP.data(:,HDR.InChanSelect);                      
                                
                        elseif HDR.SCP2.NHT~=19999,
                                fprintf(HDR.FILE.stderr,'Warning SOPEN SCP-ECG: user specified Huffman Table not supported\n');
                                return;
                        else
                                HDR.SCP2,
                        end;

                        % Decoding of Difference encoding                  
                        if SCP.FLAG.DIFF==2,
                                for k1 = 3:size(S2,1);
                                        S2(k1,:) = S2(k1,:) + [2,-1] * S2(k1-(1:2),:);
                                end;
                        elseif SCP.FLAG.DIFF==1,
                                S2 = cumsum(S2);    
                        end;
                        S2 = S2 * SCP.Cal;

			if (k3==5) && isfield(HDR,'SCP5');
			%	HDR.SCP5.data = S2; 

			elseif (k3==6) && isfield(HDR,'SCP6');
                                if HDR.SCP6.FLAG.bimodal_compression,
                                        F = HDR.SCP5.SampleRate/HDR.SCP6.SampleRate;
                                        HDR.SampleRate = HDR.SCP5.SampleRate;
                                        HDR.FLAG.F = F;
                                        
                                        tmp=[HDR.SCP4.PA(:,1);HDR.LeadPos(1,2)]-[1;HDR.SCP4.PA(:,2)+1];
                                        if ~all(tmp==floor(tmp))
                                                tmp,
                                        end;
                                        t  = (1:HDR.N) / HDR.SampleRate;
                                        S1 = zeros(HDR.N, HDR.NS);
                                        
                                        
                                        p = 1;
                                        k2 = 1;
                                        pa = [HDR.SCP4.PA;NaN,NaN];
                                        flag = 1;
                                        for k1 = 1:HDR.N,
                                                if k1 == pa(p,2)+1,
                                                        flag = 1;
                                                        p    = p+1;
                                                        accu = S2(k2,:);
                                                elseif k1 == pa(p,1),
                                                        flag = 0;
                                                        k2 = ceil(k2);
                                                end;
                                                
                                                if flag,
                                                        S1(k1,:) = ((F-1)*accu + S2(floor(k2),:)) / F;
                                                        k2 = k2 + 1/F;
                                                else	
                                                        S1(k1,:) = S2(k2,:);
                                                        k2 = k2 + 1;
                                                end;
                                        end;	
                                        
                                        HDR.SCP.S2 = S2;
                                        HDR.SCP.S1 = S1;
                                        S2 = S1;
                                end;
                                
                                if HDR.FLAG.ReferenceBeat,
                                        for k = find(~HDR.SCP4.type(:,1)'),
                                                t1 = (HDR.SCP4.type(k,2):HDR.SCP4.type(k,4));
                                                t0 = t1 - HDR.SCP4.type(k,3) + HDR.SCP4.fc0;
                                                S2(t1,:) = S2(t1,:) + HDR.SCP5.data(t0,:); 
                                        end;
                                end;
			end;
		end;
	end;
        HDR.data = S2;
	HDR.TYPE = 'native'; 	% decompression is already applied

	if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;

        nr = min(round(HDR.SampleRate * NoS), size(HDR.data,1) - HDR.FILE.POS);
        S  = HDR.data(HDR.FILE.POS + (1:nr), :);
        HDR.FILE.POS = HDR.FILE.POS + nr;
        

elseif strcmp(HDR.TYPE,'SIGIF'),
        if nargin==3,
                HDR.FILE.POS = StartPos;
        end;
        
        S = [];
        for k = 1:min(NoS,HDR.NRec-HDR.FILE.POS),
                HDR.FILE.POS = HDR.FILE.POS + 1;
                STATUS = fseek(HDR.FILE.FID, HDR.Block.Pos(HDR.FILE.POS), 'bof');
                if HDR.FLAG.TimeStamp,
                        HDR.Frame(k).TimeStamp = fread(HDR.FILE.FID,[1,9],'uint8');
                end;
                
                if HDR.FLAG.SegmentLength,
                        HDR.Block.Length(k) = fread(HDR.FILE.FID,1,'uint16');  %#26
                        STATUS = fseek(HDR.FILE.FID,HDR.Block.Length(k)*H1.Bytes_per_Sample,'cof');
                else
                        tmp = HDR.Segment_separator-1;
                        [dat,c] = fread(HDR.FILE.FID,[HDR.NS,HDR.Block.Length/HDR.NS],gdfdatatype(HDR.GDFTYP));
                        [tmpsep,c] = fread(HDR.FILE.FID,1,gdfdatatype(HDR.GDFTYP));
                        
                        if  (tmpsep~=HDR.Segment_separator);
                                fprintf(HDR.FILE.stderr,'Error SREAD Type=SIGIF: blockseparator not found\n');
                        end;
                end;
                S = [S; dat(HDR.InChanSelect,:)'];
        end;
        
        
elseif strcmp(HDR.TYPE,'CTF'),
        if nargin>2,
                STATUS = fseek(HDR.FILE.FID,HDR.HeadLen+HDR.NS*HDR.SPR*4*StartPos,'bof');        
                HDR.FILE.POS = StartPos;
        end;
	
	nr = min(NoS, HDR.NRec - HDR.FILE.POS);
	
        S = []; count = 0; 
	for k = 1:nr,
	        %[tmp,c] = fread(HDR.FILE.FID, 1, 'int32')
	        [s,c] = fread(HDR.FILE.FID, [HDR.SPR, HDR.NS], 'int32');
		S = [S; s(:,HDR.InChanSelect)];
		count = count + c;
	end;
	
        HDR.FILE.POS = HDR.FILE.POS + count/(HDR.SPR*HDR.NS);
        
        
elseif strcmp(HDR.TYPE,'EEProbe-CNT'),
        if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        
        nr = min(HDR.SampleRate*NoS, HDR.SPR*HDR.NRec-HDR.FILE.POS);
	if exist('read_eep_cnt','file')==3,
                tmp = read_eep_cnt(HDR.FileName, HDR.FILE.POS+1, HDR.FILE.POS+nr);
                sz  = size(tmp.data);
                S   = tmp.data(HDR.InChanSelect,:)';
                clear tmp; 
                if HDR.FLAG.UCAL,
	                S = S*diag(1./HDR.Cal(HDR.InChanSelect));
	        end;        
                HDR.FILE.POS = HDR.FILE.POS + sz(2);

        elseif exist('OCTAVE_VERSION','builtin')
                fprintf(HDR.FILE.stderr,'ERROR SREAD (EEProbe): Reading EEProbe-file format is not supported.\n');
                return;
	else
                fprintf(HDR.FILE.stderr,'ERROR SREAD (EEProbe): Cannot open EEProbe-file, because read_eep_cnt.mex not installed. \n');
                fprintf(HDR.FILE.stderr,'ERROR SREAD (EEProbe): You can downlad it from http://www.smi.auc.dk/~roberto/eeprobe/\n');
		%% ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/external/eeprobe.zip
                return;
        end
       
elseif strcmp(HDR.TYPE,'EEProbe-AVR'),
        if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;

        nr = min(HDR.SPR-HDR.FILE.POS,NoS*HDR.SampleRate);

        S = HDR.EEP.data(HDR.FILE.POS+(1:nr),:);
        HDR.FILE.POS = HDR.FILE.POS + nr;

        
elseif strncmp(HDR.TYPE,'BrainVision',11),   %Brainvision
        if strncmpi(HDR.BV.DataFormat, 'binary',5)
		tc = strcmp(HDR.TYPE,'BrainVisionVAmp');
		NS = HDR.NS+tc;
                if strncmpi(HDR.BV.DataOrientation, 'multiplexed',6),
                        if nargin>2,
                                STATUS = fseek(HDR.FILE.FID,StartPos*HDR.SampleRate*HDR.AS.bpb,'bof');        
                                HDR.FILE.POS = HDR.SampleRate*StartPos;
                        end;

			nr = min(HDR.SampleRate*NoS, HDR.SPR*HDR.NRec - HDR.FILE.POS);
			if (length(HDR.InChanSelect)*2>HDR.NS)
				[s,c] = fread(HDR.FILE.FID, [NS, nr], ['*',gdfdatatype(HDR.GDFTYP)]);
				count = c/NS;
				S = s(HDR.InChanSelect,:)';
			else
				S  = [];
				count = 0;
				while (count<nr),
					[s,c] = fread(HDR.FILE.FID, [NS, min(nr-count,floor(2^24/NS))], gdfdatatype(HDR.GDFTYP));
					if ~c, break; end; 
					S = [S; s(HDR.InChanSelect,:)'];
					count = count + c/NS; 
				end; 
			end;
			HDR.FILE.POS = HDR.FILE.POS + count;
		elseif strncmpi(HDR.BV.DataOrientation, 'vectorized',6),
			S = [];
			nr = min(HDR.SampleRate*NoS, HDR.AS.endpos-HDR.FILE.POS);

			count = 0; 
			for chan = 1:length(HDR.InChanSelect);
                                STATUS = fseek(HDR.FILE.FID, HDR.HeadLen + HDR.FILE.POS + HDR.AS.bpb*HDR.SPR*(chan-1)/NS, 'bof');
                                [s,count] = fread(HDR.FILE.FID, [nr,1], gdfdatatype(HDR.GDFTYP));
                                if count ~= nr,
                                        fprintf(2,'ERROR READ BV-bin-vec: \n');
                                        return;
                                end;
                                S(:,chan) = s;
                        end
                        HDR.FILE.POS = HDR.FILE.POS + count; 
                end;

        elseif strncmpi(HDR.BV.DataFormat, 'ascii',5)  
	        %%%% OBSOLETE: supported by 'native' %%%%
                if nargin>2,
                        HDR.FILE.POS = HDR.SampleRate*StartPos;
                end;
                nr = min(HDR.SampleRate*NoS, HDR.AS.endpos-HDR.FILE.POS);
                S  = HDR.BV.data(HDR.FILE.POS+(1:nr),HDR.InChanSelect);
                
        end
                
        
elseif strcmp(HDR.TYPE,'SierraECG'),   %% SierraECG  1.03  *.open.xml from PHILIPS
        if ~isfield(HDR,'data');
                [HDR.data,status] = str2double(HDR.XML.waveforms.parsedwaveforms);
                if any(status)
                        error('SREAD: compressed SierraECG (Philips) format not supported')
                end;
                HDR.data = reshape(HDR.data,length(HDR.data)/HDR.NS,HDR.NS);
                HDR.SPR = size(HDR.data,1);
        else
                % base64 - decoding 
                base64 = ['A':'Z','a':'z','0':'9','+','/'];
                decode64 = repmat(nan,256,1);
                decode64(abs(base64)) = 0:63;
                tmp = decode64(HDR.XML.waveforms.parsedwaveforms);
                tmp(isnan(tmp)) = [];
                n   = length(tmp);
                tmp = reshape([tmp;zeros(mod(n,4),1)], 4, ceil(n/4));
                t1  = tmp(1,:)*4 + floor(tmp(2,:)/16);
                t2  = mod(tmp(2,:),16)*16 + floor(tmp(3,:)/4);
                t3  = mod(tmp(3,:),4)*64 + tmp(4,:);
                tmp = reshape([t1,t2,t3], ceil(n/4)*3, 1);
        end;
        if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        nr = min(HDR.SampleRate*NoS, HDR.SPR-HDR.FILE.POS);
        S  = HDR.data(HDR.FILE.POS+(1:nr),HDR.InChanSelect);
        HDR.FILE.POS = HDR.FILE.POS + nr;

        
elseif strcmp(HDR.TYPE,'ATF'); 
        if HDR.FILE.OPEN,
                fseek(HDR.FILE.FID,HDR.HeadLen,-1);
                t = fread(HDR.FILE.FID,[1,inf],'uint8');
                fclose(HDR.FILE.FID);
                HDR.FILE.OPEN=0; 
                [HDR.ATF.NUM,status,HDR.ATF.STR] = str2double(char(t));
        end;
        S = HDR.ATF.NUM;
        
        
elseif strcmp(HDR.TYPE,'FEPI3'); 
        if nargin==3,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        Duration = min(NoS*HDR.SampleRate,(HDR.AS.endpos-HDR.FILE.POS));

	S = repmat(NaN,Duration,length(HDR.InChanSelect));
	ix = find([HDR.FEPI.SEG(:,2) > HDR.FILE.POS] & [HDR.FEPI.SEG(:,1) <= HDR.FILE.POS+Duration]);
	for k = 1:length(ix),
		fid = fopen(fullfile(HDR.FILE.Path,[HDR.FEPI.ListOfDataFiles{k},'.bin']));
		if k==1,
			fseek(fid,(HDR.FILE.POS-HDR.FEPI.SEG(ix(k),1))*2*HDR.NS,-1);
		end;	
		data= fread(fid,[HDR.NS,inf],'int16')'; 
		ix1 = HDR.FEPI.SEG(ix(k),1)-HDR.FEPI.SEG(ix(1),1);
		ix2 = min(size(data,1),Duration+HDR.FILE.POS-HDR.FEPI.SEG(ix(k),1)+1);
	        S(ix1+1:ix1+ix2,:) = data(1:ix2,HDR.InChanSelect);
		fclose(fid); 
	end; 
	HDR.FILE.POS = HDR.FILE.POS+Duration; 
        
        
elseif strcmp(HDR.TYPE,'WG1'),   %walter-graphtek
	% code from Robert Reijntjes, Amsterdam, NL 
	% modified by Alois Schloegl 19. Feb 2005 
        if nargin==3,
                HDR.FILE.POS = round(HDR.SampleRate*StartPos);
        end;

	ix1    = mod(HDR.FILE.POS, HDR.SPR);	% starting sample (minus one) within 1st block 
    	fp     = HDR.HeadLen + floor(HDR.FILE.POS/HDR.SPR)*HDR.AS.bpb;
    	status = fseek(HDR.FILE.FID, fp, 'bof');

        nr     = min(HDR.AS.endpos-HDR.FILE.POS, NoS*HDR.SampleRate);
	S      = repmat(NaN,nr,length(HDR.InChanSelect)); 
	count  = 0;
        endloop= 0;
        c = 1; 
        offset = 0; 
    	while ~endloop & (c>0) && (offset(1)~=(hex2dec('AEAE5555')-2^32)) && (count<nr);
	        [offset,c] = fread(HDR.FILE.FID, HDR.WG1.szOffset, 'int32');
		[databuf,c] = fread(HDR.FILE.FID,[HDR.WG1.szBlock,HDR.NS+HDR.WG1.szExtra],'uint8');
            	dt = HDR.WG1.conv(databuf(:,1:HDR.NS)+1);
            	if any(dt(:)==HDR.WG1.unknownNr),
            		%dt(dt==HDR.WG1.unknownNr) = NaN; 
            		dt(:) = NaN; 
            	    	%fprintf(HDR.FILE.stderr,'Warning SREAD (WG1): error in reading datastream');
            	end;
		dt(1,:) = dt(1,:) + offset(1:HDR.NS)';
		dt = cumsum(dt,1);
		
		ix2 = min(nr-count, size(dt,1)-ix1);
		S(count+1:count+ix2,:) = dt(ix1+1:ix1+ix2, HDR.InChanSelect);
		count = count + ix2; 
		ix1 = 0;	% reset starting index, 

                k = 0; 
                while (k<HDR.WG1.szExtra) && ~endloop, 
                        endloop = ~isempty(strfind(databuf(:,HDR.NS+k)',[85,85,174,174]));
                        k = k+1; 
                end;
	end;	
	%S = S(1:count,:);
	HDR.FILE.POS = HDR.FILE.POS + count;


elseif strcmp(HDR.TYPE,'XML-FDA'),   % FDA-XML Format
        if ~isfield(HDR,'data');
                tmp   = HDR.XML.component.series.derivation;
                if isfield(tmp,'Series');
                        tmp = tmp.Series.component.sequenceSet.component;
                else    % Dovermed.CO.IL version of format
                        tmp = tmp.derivedSeries.component.sequenceSet.component;
                end;
                for k = 1:length(HDR.InChanSelect);
                        HDR.data(:,k) = str2double(tmp{HDR.InChanSelect(k)+1}.sequence.value.digits)';
                end;
                HDR.SPR = size(HDR.data,1);
        end;
        if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;
        nr = min(HDR.SampleRate*NoS, HDR.SPR-HDR.FILE.POS);
        S  = HDR.data(HDR.FILE.POS+(1:nr),:);
        HDR.FILE.POS = HDR.FILE.POS + nr;

% using XML4MAT instead of XMLTREE
% str2double(HDR.XML0{end}.component{1}.series{end}.derivation{:}.derivedSeries{5}.component{1}.sequenceSet{5}.component{1}.sequence{2}.value{3}.digits)'

        
elseif strcmp(HDR.TYPE,'FIF'),
        % some parts of this code are from Robert Oostenveld, 
        if ~(exist('rawdata')==3 & exist('channames')==3)
                error('cannot find Neuromag import routines on your Matlab path (see http://boojum.hut.fi/~kuutela/meg-pd)');
        end
        if nargin<3, 
                StartPos = HDR.FILE.POS/HDR.SampleRate;
        end
        if nargin>2,
                HDR.FILE.POS = HDR.SampleRate*StartPos;
        end;

        t1  = rawdata('goto', HDR.FILE.POS/HDR.SPR);
        t2  = t1;
        dat = [];
        count = 0;
        status = 'ok';
        
        while (t2<(StartPos + NoS)) && ~strcmp(status,'eof'),
                [buf, status] = rawdata('next');
                if 0
                elseif strcmp(status, 'ok')
                        count = count + size(buf,2);
                        dat = [dat; buf(HDR.InChanSelect,:)'];
                elseif strcmp(status, 'eof')
                elseif strcmp(status, 'skip')
                elseif strcmp(status, 'error')
                        error('error reading selected data from fif-file');
                else
                        error('undefined status code return from RAWDATA(FIF-file)');
                end
                t2 = rawdata('t');
        end
        t  = t1*HDR.SampleRate+1:t2*HDR.SampleRate;
        ix = (t>StartPos*HDR.SampleRate) & (t<=(StartPos+NoS)*HDR.SampleRate);
        S  = dat(ix,HDR.InChanSelect);
        HDR.FILE.POS = t2*HDR.SampleRate;        

elseif strcmp(HDR.TYPE,'EVENT'),
        s = [];        

elseif strncmp(HDR.TYPE,'IMAGE:',6),
	% forward call to IREAD
        [S,HDR] = iread(HDR);
	return;

else
        fprintf(2,'Error SREAD: %s-format not supported yet.\n', HDR.TYPE);        
	return;
end;


%%% TOGGLE CHECK - checks whether HDR is kept consist %%% 
if 0,
global SREAD_TOGGLE_CHECK
if isfield(HDR.FLAG,'TOGGLE');
        if HDR.FLAG.TOGGLE~=SREAD_TOGGLE_CHECK,
                fprintf(HDR.FILE.stderr,'Warning SREAD: [s,HDR]=sread(HDR, ...) \nYou forgot to pass HDR in %i call(s) of SREAD\n',SREAD_TOGGLE_CHECK-HDR.FLAG.TOGGLE);
        end;
else
        HDR.FLAG.TOGGLE=0;
        SREAD_TOGGLE_CHECK=0;
end;
SREAD_TOGGLE_CHECK = SREAD_TOGGLE_CHECK+1;
HDR.FLAG.TOGGLE = HDR.FLAG.TOGGLE+1;
end; 

if STATUS,
        fprintf(HDR.FILE.stderr,'WARNING SREAD: something went wrong. Please send the files %s and BIOSIGCORE to <a.schloegl@ieee.org>',HDR.FileName);
        save biosigcore.mat 
end;

if isempty(S),

elseif isfield(HDR,'THRESHOLD') && HDR.FLAG.OVERFLOWDETECTION,
        ix = (S~=S);
        for k=1:length(HDR.InChanSelect),
                TH = THRESHOLD(HDR.InChanSelect(k),:);
                %ix(:,k) = (S(:,k)<=TH(1)) | (S(:,k)>=TH(2));
                ix = (S(:,k)<=TH(1)) | (S(:,k)>=TH(2));
                S(ix,k)=NaN;
        end
        if exist('double','builtin')
                S = double(S);
        end;
%        S(ix>0) = NaN;
elseif HDR.FLAG.OVERFLOWDETECTION,
        % no HDR.THRESHOLD defined
        warning('no Threshold defined');
elseif isfield(HDR,'THRESHOLD'),
        % automated overflow detection has been turned off
end;

if ~HDR.FLAG.UCAL,
        % S = [ones(size(S,1),1),S]*HDR.Calib; 
        % perform the previous function more efficiently and
        % taking into account some specialities related to Octave sparse
        % data. 
        if isempty(S),	% otherwise Octave 2.1.64 could break below, 
		if size(S,2)~=length(HDR.InChanSelect), 
			fprintf(HDR.FILE.stderr,'Warning SREAD (%s): number of columns (%i) incorrect!\n',HDR.TYPE,size(S,2));
		end;	
		S = zeros(0,length(HDR.InChanSelect));
	end;

        %if ~issparse(HDR.Calib); %
        if FLAG_CALIB_DONE, 

        elseif strcmpi(HDR.FLAG.OUTPUT,'single')
        	tmp = single(zeros(size(S,1),size(HDR.Calib,2)));
                for k = 1:size(S,1),
                        tmp(k,:) = [1,S(k,:)] * HDR.Calib;
                end;
                S = tmp; 
                clear tmp;
        	
        elseif 1, % exist('OCTAVE_VERSION','builtin')
                % force octave to do a sparse multiplication
                % the difference is NaN*sparse(0) = 0 instead of NaN
                % this is important for the automatic overflow detection
		Calib = HDR.Calib;
		tmp   = S;
                S     = zeros(size(S,1),size(Calib,2));   % memory allocation

                for k = 1:size(Calib,2),
                        chan = find(Calib(2:end,k));
                        S(:,k) = double(tmp(:,chan)) * full(Calib(1+chan,k)) + Calib(1,k);
                end;
        else
                % S = [ones(size(S,1),1),S]*HDR.Calib; 
                % the following is the same as above but needs less memory. 
                S = full(double(S) * HDR.Calib(2:end,:));
                for k = 1:size(HDR.Calib,2),
                        S(:,k) = S(:,k) + full(HDR.Calib(1,k));
                end;
        end;
end;



