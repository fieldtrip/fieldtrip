function [HDR,data,t]=matread(HDR,arg2,idxlist)
% MATRREAD Loads (parts of) data stored in Matlab-format 
%
% [HDR,data,timeindex]=matread(HDR,block_number, [startidx, endidx])
% This is the recommended use for Matlab-files generated from ADICHT data
% Before using MATREAD, HDR=MATOPEN(filename, 'ADI', ...) must be applied.
%
% [HDR,data,timeindex]=matread(HDR,Variable_Name, [startidx, endidx])
% can be used for other Matlab4 files. 
% Variable name is a string which identifies a Matlab Variable.
% Before using MATREAD, HDR=MATOPEN(filename, 'r', ...) must be applied.
%
% see also: EEGREAD, FREAD, EEGOPEN, EEGCLOSE  

%	$Revision$
%	$Id$
%	Copyright (c) 1997-2003 by  Alois Schloegl
%	a.schloegl@ieee.org	

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

CHt=[];
if ~ischar(arg2)
        if strcmp(HDR.TYPE,'ADI'); %HDR.ADI.Mode,
                BlockNo=abs(arg2);
                if isempty(BlockNo), BlockNo=1:length(HDR.ADI.DB); end;
                if (isinf(BlockNo) | (length(BlockNo)~=1))
                        if nargout<3,
		                fprintf(2,'Warning MATREAD: missing output argument, time information cannot be returned\n');
                        end;
                        CHd=HDR.ADI.DB(BlockNo);
                        CHt=HDR.ADI.TB(BlockNo);
                elseif ((BlockNo<1) | BlockNo>length(HDR.ADI.DB)),
	                fprintf(2,'Warning MATREAD: tried to read block %i, but only %i blocks in file %s\n',BlockNo,length(HDR.ADI.DB),HDR.FileName);
                        CHd=[];CHt=[];
                else
                        CHd=HDR.ADI.DB(BlockNo);
                        CHt=HDR.ADI.TB(BlockNo);
                end;
        else 
                CHd=arg2;
                fprintf(2,'Warning MATREAD: Variable identified by number\n');
        end;	
else  %if ischar(arg2)
        CHd=find(strcmp(arg2,{HDR.Var.Name}));
        if isempty(CHd)
                %HDR.ErrNo=-1; 
                fprintf(2,'Warning MATREAD: Variable %s not found\n',arg2);
                data=[];
                return;
        end;
end;
if length(CHd)>1 
        fprintf(2,'Warning MATREAD: Only one block/variable can be read\n');
        CHd=CHd(1);
end;

if nargin<4, end;

if nargin<3,
        idxlist=[1,HDR.Var(CHd).Size(2)];	        
else
	if ~all(isfinite(idxlist)),
        	idxlist=[1,HDR.Var(CHd).Size(2)];	        
	end;        
end

% if length(idxlist)<2,idxlist=idxlist*[1 1]; end;

if (idxlist(1)<1) %| (idxlist(length(idxlist))>HDR.Var(CHd).Size(2))), 
        fprintf(2,'Warning MATREAD #%i: Invalid File Position %f-%f\n',CHd,[idxlist(1)-1,idxlist(length(idxlist))]); 
        idxlist(1)=1;
        %return; 
end;
if idxlist(length(idxlist))>HDR.Var(CHd).Size(2), 
        fprintf(2,'Warning MATREAD #%i: endidx exceeds block length %f-%f\n',CHd,[HDR.Var(CHd).Size(2),idxlist(length(idxlist))]); 
        idxlist=[idxlist(1),HDR.Var(CHd).Size(2)];
        %return; 
end;
if 0; %((idxlist(1)<1) | (idxlist(length(idxlist))>HDR.Var(CHd).Size(2))), 
        fprintf(2,'ERROR MATREAD #%i: Invalid File Position %f-%f\n',CHd,[idxlist(1)-1,idxlist(length(idxlist))]); 
        return; 
end;
Pos = (idxlist(1)-1) * HDR.Var(CHd).Size(1) * HDR.Var(CHd).SizeOfType;
Len = min(idxlist(length(idxlist)), HDR.Var(CHd).Size(2)) - idxlist(1) + 1;

fseek(HDR.FILE.FID, round(HDR.Var(CHd).Pos + Pos), -1); % round must be explicite, otherwise fseek does incorrect rounding
[msg,errno] = ferror(HDR.FILE.FID); 
if errno, fprintf(2,'ERROR MATREAD: FSEEK does not work Pos=%f\n%s',HDR.Var(CHd).Pos+Pos,msg); return; end;

count = 0;

if ~any(CHd==HDR.ADI.DB), % for non-data blocks 
        data=repmat(nan,HDR.Var(CHd).Size(1),Len);
	while (count<Len),
        	cc = min(2^12,Len-count);
		% [CHd,idxlist(1),idxlist(length(idxlist)),Len,HDR.Var(CHd).Size(2)]
		[dta,c] = readnextblock4(HDR.FILE.FID,HDR.Var(CHd),cc);
                data(:,count+(1:size(dta,2))) = dta;
	        count = count + cc;
	end;
else    
        BLOCKSIZE=64*3*25*4; % 2^12; % You can change this for optimizing on your platform
        
        %HDR=mat_setfilter(HDR,BlockNo);% set filters for this block 
        for k=1:max(HDR.SIE.InChanSelect),
                if k <= HDR.NS,
		        if HDR.SIE.FILT==1;
                                [tmp,HDR.Filter.Z(:,k)] = filter(HDR.Filter.B,HDR.Filter.A,zeros(length(HDR.Filter.B),1));
		                HDR.FilterOVG.Z = HDR.Filter.Z;
                        end;
                end;
        end;
        
        if isfield(HDR,'iFs') & (HDR.iFs>0) & (HDR.iFs<inf),                        
                Fs=HDR.SampleRate(find(CHd==HDR.ADI.DB))/HDR.iFs;
	        data = repmat(nan,ceil(Len*size(HDR.SIE.T,2)/size(HDR.SIE.T,1)/Fs),length(HDR.SIE.ChanSelect));
        else
                Fs=1;
	        data = repmat(nan,ceil(Len*size(HDR.SIE.T,2)/size(HDR.SIE.T,1)),length(HDR.SIE.ChanSelect));
        end;
	count=0; count2=0;
        while (count2<Len), 
                cc = min(BLOCKSIZE,Len-count2);
                % [Len, count2, count, Fs, (Len-count2)*Fs, cc]
                % [CHd,idxlist(1),idxlist(length(idxlist)),Len,HDR.Var(CHd).Size(2)]
                [dta,c] = readnextblock4(HDR.FILE.FID,HDR.Var(CHd),cc);
                % ferror(HDR.FILE.FID)
                % [Len,cc,count,BLOCKSIZE,HDR.FILE.Pos,ftell(HDR.FILE.FID)],
                [nc,cc] = size(dta);
                
                % if any postprocessing, resample to internal Sampling rate iFs=1000Hz)
                if Fs~=1,
                        dta=rs(dta',Fs,1)'; %**************************************************************************************1
                end;

                if nc~=HDR.NS, % reorganizing of channels
	                dta=sparse(find(HDR.ADI.index{arg2}),1:sum(HDR.ADI.index{arg2}>0),1)*dta;
        	        dta(find(~HDR.ADI.index{arg2}),:)=nan;
                end;
                
                for k = sort(HDR.SIE.InChanSelect), %1:size(data,2);
                                %k=HDR.SIE.chanselect(K),
                	if k<=HDR.NS,
		                if HDR.SIE.FILT,
                                        [dta(k,:),HDR.Filter.Z(:,k)]=filter(HDR.Filter.B,HDR.Filter.A,dta(k,:),HDR.Filter.Z(:,k));
				end;                                        
                        end;
                end;	
                
                if HDR.SIE.RS,
                        dta = rs(dta(HDR.SIE.ChanSelect,:)',HDR.SIE.T); %RS%              ***********************************************2
                else
                    	dta = dta(HDR.SIE.ChanSelect,:)';    
                end;
                data(count+(1:size(dta,1)),:) = dta;
	        count  = count  + size(dta,1);
	        count2 = count2 + cc;
	end;
        
        %%%%% Overflow Detection %%%%%
        if HDR.SIE.TH,
                for k=1:find(~isnan(sum(HDR.SIE.THRESHOLD,2))),
	                ch  = phrchan(HDR.FILE.Name);
        	        tmp = (data(:,k)<HDR.SIE.THRESHOLD(1,k)) | (data(:,ch)>HDR.SIE.THRESHOLD(2,k));
                        data(tmp,k) = NaN;
                end;	
        end;
        
        % load time vector
	if nargout>2,
	if ~isempty(CHt),
		fseek(HDR.FILE.FID, round(HDR.Var(CHt).Pos + (idxlist(1)-1)*HDR.Var(CHt).SizeOfType),-1);
                [t,c]=readnextblock4(HDR.FILE.FID,HDR.Var(CHt),idxlist(length(idxlist))-min(idxlist)+1);
                if length(t)>1,
	                tmp=diff(t); 
        	        if any(abs(tmp-tmp(length(tmp)))>1000*eps)
                                fprintf(2,'Warning %s: Sampling is not equally spaced in %s [%i,%i]\n',mfilename,HDR.FileName,idxlist(1),idxlist(length(idxlist)));
                        end;	
                end;
                
		%%%%% Delay of Filtering
                if HDR.SIE.FILT,
                        t = t' + HDR.Filter.Delay;
                else
                        t = t';
                end;
                
                % if any postprocessing, resample to internal Sampling rate iFs=1000Hz)
                if isfield(HDR,'iFs') & (HDR.iFs>0) & (HDR.iFs<inf),                        
	                Fs=HDR.SampleRate(find(CHd==HDR.ADI.DB))/HDR.iFs;
        	        if Fs~=1,
                	        t=rs(t,Fs,1);           %*********************************************************************3
	                end;     
                end;
                
                %%%%% Resampling of the time
		if HDR.SIE.RS, 
                        t=rs(t,HDR.SIE.T); %RS%    %******************************************************************4
                end;
	else
       	       	fprintf(2,'Warning %s: timeindex is not available \n',mfilename);
	end;

	end;
end;


function [data,c]=readnextblock4(fid,VarInfo,Len),
       	dt=VarInfo.Type(3);
       	if     dt==0,       type = 'float64';
       	elseif dt==6,       type = 'uint8';
       	elseif dt==4,       type = 'uint16';
       	elseif dt==3,       type = 'int16';
       	elseif dt==2,       type = 'int32';
       	elseif dt==1,       type = 'float32';
        else
               fprintf(2,'Error %s: unknown data type\n',mfilename);
               return;
	end;
                
        [data,c]=fread(fid,[VarInfo.Size(1),Len],type);

        if VarInfo.Type(5);
                %HDR.ErrNo=-1;
                fprintf(2,'Warning %s: imaginary data not test\n',mfilename);
                        
                [di,c]=fread(fid,[VarInfo.Size(1),Len],type);
                data=data+i*di;
      	end;

	if VarInfo.Type(4)==1,data=char(data); end;
        
        