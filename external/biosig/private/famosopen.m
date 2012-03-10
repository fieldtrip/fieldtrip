function [HDR]=famosopen(arg1,arg3,arg4,arg5,arg6)
% FAMOSOPEN opens FAMOS file
% However, it is recommended to use SOPEN instead .
% For loading whole data files, use SLOAD. 
%
% see also: SOPEN, SREAD, SSEEK, STELL, SCLOSE, SWRITE, SEOF

% HDR=famosopen(HDR);

%	$Id$
%	Copyright (c) 2007 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if isstruct(arg1),
	HDR=arg1;
else
	FILENAME=arg1;
        HDR.FileName  = FILENAME;
        [pfad,file,FileExt] = fileparts(HDR.FileName);
        HDR.FILE.Name = file;
        HDR.FILE.Path = pfad;
        HDR.FILE.Ext  = FileExt(2:length(FileExt));
	HDR.FILE.PERMISSION = 'r';
	HDR.FILE.stdout = 1; 
	HDR.FILE.stderr = 2; 
end;


if any(HDR.FILE.PERMISSION=='r'),
        %%%%% READ HEADER
        HDR.FILE.FID = fopen(HDR.FileName,'r');
        if HDR.FILE.FID<0,
                fprintf('Error SOPEN: could not open file %s\n',HDR.FileName);
                return;
        end; 
        POS = 0; 
        tline = [];
        CHAN = 0; 
        HDR.AS.SampleRate =[];
        while ~feof(HDR.FILE.FID), 
        	%POS = ftell(HDR.FILE.FID);
	        [t,count] = fread(HDR.FILE.FID,[1,1024],'*char');
	        tline = [tline,t];
		
		while length(tline)>512
			% extract key
			ix = find(tline=='|'); ix1=ix(1); 
			tline = tline(ix1+1:end); 
		        POS = POS + ix1;
			
			% extract key, version, length
			ix = find(tline==','); ix = ix(1:3); 

		        [n,v,s] = str2double(tline(1:ix(3)),',');
		        keylen  = n(3); 
		        if ix(3)+1+keylen <= length(tline),
			        param  = tline(ix(3)+1:ix(3)+1+keylen);
			else 
				param = tline(ix(1)+1:35);        
			end; 
		        tline = tline(ix(3)+1:end);
		        POS = POS + ix(3);
		
	        	fprintf(1,'SOPEN(FAMOS) %i: process "%s %s"\n',POS,s{1},param);
		        if isempty(s),
		        elseif strcmp(s{1},'CF') & (n(2)==2)
		        	if ~all(n(2:3)==[1,1])
		        		fprintf(HDR.FILE.stdout,'Warning SOPEN(FAMOS): unknown/unsupported(?) version\n');
		        	end; 	
		        	
		        elseif strcmp(s{1},'CK') & (n(2)==1)
		        	if ~all(n(2:3)==[3,1])
		        		fprintf(HDR.FILE.stdout,'Warning SOPEN(FAMOS): file %s was not closed, it is perhaps corrupted\n',HDR.FILE.Name);
		        	end; 	
		        	
		        elseif 0, strcmp(s{1},'NO') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif 0, strcmp(s{1},'CT') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif 0, strcmp(s{1},'CB') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif 0, strcmp(s{1},'CI') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif 0, strcmp(s{1},'CG') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif 0, strcmp(s{1},'CD') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif 0, strcmp(s{1},'CC') & (n(2)==1)
			        [n1,v1,s1] = str2double(param,',');
%			        CHAN = CHAN+1;
		        	
		        elseif 0, strcmp(s{1},'NT') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif 0, strcmp(s{1},'CZ') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif 0, strcmp(s{1},'CC') & (n(2)==1)
			        %[n1,v1,s1] = str2double(param,',')
		        	
		        elseif strcmp(s{1},'CP') & (n(2)==1)
			        [n1,v1,s1] = str2double(param,',');
			        CHAN = n1(1); 
				switch n1(3),
				case 1, typ=2;	% uint8
				case 2, typ=1;  % int8
				case {3,9,11}, typ=4;	% uint16
				case 4, typ=3;  % int16
				case 5, typ=6;	% uint32
				case 6, typ=5;	% int32   
				case 7, typ=16;	% float
				case 8, typ=17;	% double
				case 10,typ=0;	% double
				case 13,typ=511+48;	% double
				otherwise,
			                fprintf(1,'SOPEN(FAMOS): typ not supported (yet)\n');
			                return;
				end; 				 	
				HDR.GDFTYP(CHAN)=typ;
				if n1(5),
			                fprintf(1,'SOPEN(FAMOS): mask %i not supported!\n',n1(5));
				end; 			        
				if n1(5:8)~=[0,0,1,0],
			                fprintf(1,'SOPEN(FAMOS): mask %i.%i.%i.%i not supported (yet)\n',n1(5:8));
				end; 			        
		        	
		        elseif strcmp(s{1},'Cb') & (n(2)==1)
			        [n1,v1,s1] = str2double(param,',');
			        ch = n1(3); 
				sz = n1(2); % bytes per sample
				HDR.AS.start(ch) = n1(5);
				HDR.AS.bytes(ch) = n1(6);
		        	
		        elseif 0, strcmp(s{1},'CG') & (n(2)==1)
			        [n1,v1,s1] = str2double(param,',');
		        	
		        elseif strcmp(s{1},'CD') & (n(2)==1)
			        [n1,v1,s1] = str2double(param,',');
			        HDR.AS.SampleRate(length(HDR.AS.SampleRate)+1) = 1/n1(1);
			        HDR.AS.PhysDimDeltaT = s1{4};
		        	
		        elseif strcmp(s{1},'CR') & (n(2)==1)
			        [n1,v1,s1] = str2double(param,',');
				if length(s1)<6,			        
				        HDR.PhysDim{CHAN} = ' '; 
			        else
				        HDR.PhysDim{CHAN} = s1{6}; 
				end;        
			        HDR.LeadIdCode(CHAN) = 0;
			        %HDR.Cal(CHAN) = 1; 
		        	
		        elseif strcmp(s{1},'CN') & (n(2)==1)
			        [n1,v1,s1] = str2double(param,',')
			        HDR.Label{CHAN} = s1{5}; 
		        	
		        elseif strcmp(s{1},'CS') & (n(2)==1)
				ix1 = find(tline==','); ix1=min(ix1);
		        	HDR.FAMOS.LEN = n(3)- ix1; 
		        	HDR.FAMOS.POS = POS + ix1; 

		        	% read data
		        	for ch = 1:length(HDR.Label),
			        	fseek(HDR.FILE.FID,HDR.FAMOS.POS+HDR.AS.start(ch),'bof');
			        	[datatyp,limits,datatypes,numbits,GDFTYP]=gdfdatatype(HDR.GDFTYP(ch));
					[d,c] = fread(HDR.FILE.FID,[HDR.AS.bytes(ch)*8/numbits,1],datatyp);
					
					HDR.data{ch} = d;
				end;
			        
		        	POS = HDR.FAMOS.POS+HDR.FAMOS.LEN+1;
		        	fseek(HDR.FILE.FID,POS,'bof');
			        [tline,count] = fread(HDR.FILE.FID,[1,1024],'*char');
		        	
		        	%idx = s(4); 
		        else 
		                fprintf(1,'SOPEN(FAMOS) %i: key "%s %s" not supported (yet)\n',POS,s{1},param);
		        end;

%		        tline = tline2;
		end;        
	end; 
	HDR.NS = CHAN;
	HDR.FILE.OPEN=1; 
	HDR.FILE.POS =0;
end; 
