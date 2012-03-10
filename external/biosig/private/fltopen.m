function [HDR]=fltopen(arg1,arg3,arg4,arg5,arg6)
% FLTOPEN opens FLT file
% However, it is recommended to use SOPEN instead .
% For loading whole data files, use SLOAD. 
%
% see also: SOPEN, SREAD, SSEEK, STELL, SCLOSE, SWRITE, SEOF

% HDR=fltopen(HDR);

%	$Id$
%	Copyright (c) 2006,2007,2008,2009 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.

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

	fid = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.hdr']),'rt');
	[r,c] = fread(fid,[1,inf],'uint8'); 
	fclose(fid); 

	r = char(r); 
	HDR.H1 = r;
	HDR.SampleRate = 1;
	FLAG.BioSig = 0; FLAG.LockBioSig = 0;
	%HDR.FLT.remark = 'Do Not Modify !!!';
	HDR.FLT = []; 

	r0 = r; 
	hdr = '';
	while any(r0=='['); %~strcmp(hdr,'System'),
		[t1,r0] = strtok(r0,'['); 	% search begin of section
		[hdr,body] = strtok(t1,']'); 	% search end of section header
		body = strtok(body,'[]');
		b = [];

		while any(body=='=')
			[tline,body] = strtok(body,[10,13]); 
			[tok1,left1] = strtok(tline,['=']); 
			[tok2,left2] = strtok(left1,[10,13,'=']);
			if 0, 
			elseif strncmp(tok1,'*',1);
			elseif ~any(tline=='='),
			elseif any(tok2=='{');
				[tok,block]=strtok(body,'{}');
				b = setfield(b,tok1,tok);
				body = body(min(find(body=='}'))+1:end);
			else
				[num,v,sa] = str2double(tok2);
				if v,
					b = setfield(b,tok1,tok2);
				else
					b = setfield(b,tok1,num);
				end; 	
			end; 
		end;
		HDR.FLT = setfield(HDR.FLT,hdr,b);
	end; 

	HDR.SPR = HDR.FLT.Dataformat.number_of_samples; 
	HDR.NRec = 1; 
	HDR.NS = HDR.FLT.System.number_of_channels; 
	[n,v,sa]=str2double(HDR.FLT.System.parameter_of_sensors); 
	%HDR.FLT.System = rmfield(HDR.FLT.System,'parameter_of_sensors');
	HDR.FLT.sensors.id = n(:,1); 
	HDR.FLT.sensors.name = sa(:,2); 
	HDR.FLT.sensors.type = n(:,3); 
	HDR.FLT.sensors.mod  = n(:,4); 
	HDR.FLT.sensors.XYZabcArea = n(:,5:11); 

%	HDR.FLT.sensors.id(n(:,1)+1,1) = n(:,1); 
%	HDR.FLT.sensors.name(n(:,1)+1,:) = sa(:,2); 
%	HDR.FLT.sensors.type(n(:,1)+1,:) = n(:,3); 
%	HDR.FLT.sensors.mod(n(:,1)+1,:)  = n(:,4); 
%	HDR.FLT.sensors.XYZabcArea(n(:,1)+1,:) = n(:,5:11); 

	[n,v,sa]=str2double(HDR.FLT.System.parameter_of_groups);
%	HDR.FLT.System = rmfield(HDR.FLT.System,'parameter_of_groups');
	HDR.FLT.groups.id = n(:,1); 
	HDR.FLT.groups.usage = n(:,2); 
	HDR.FLT.groups.name = sa(:,3); 
	HDR.FLT.groups.unit = sa(:,4); 
	HDR.FLT.groups.unit_exp = n(:,5); 
	HDR.FLT.groups.calib = n(:,6); 
	Cal_Group(n(:,1)+1)  = n(:,6).*10.^n(:,5);
	PhysDim_Group(n(:,1)+1) = sa(:,4);

	[n,v,sa]=str2double(HDR.FLT.System.parameter_of_modules); 
	HDR.FLT.System = rmfield(HDR.FLT.System,'parameter_of_modules');
	HDR.FLT.modules.id = n(:,1); 
	HDR.FLT.modules.name = sa(:,2); 
	HDR.FLT.modules.XYZabc = n(:,3:8); 
	HDR.FLT.modules.unit = n(:,9); 
	HDR.FLT.modules.exp = n(:,10); 
	HDR.FLT.modules.unitname = sa(:,11); 

	[tline,tch] = strtok(HDR.FLT.System.parameter_of_channels,[10,13]); 
	HDR.FLT.System = rmfield(HDR.FLT.System,'parameter_of_channels');
	K = 0; 
	HDR.FLT.channels.num = repmat(NaN,HDR.NS,9); 
	HDR.FLT.channels.seq = repmat(NaN,HDR.NS,1); 
	HDR.FLT.channels.id = repmat(NaN,HDR.NS,1); 
	HDR.FLT.channels.usage = repmat(NaN,HDR.NS,1); 
	HDR.FLT.channels.cal = repmat(NaN,HDR.NS,1); 
	HDR.FLT.channels.Cal = sparse(HDR.NS,HDR.FLT.System.number_of_sensors); 
	HDR.FLT.channels.grd_mode = repmat(NaN,HDR.NS,1); 
	HDR.FLT.channels.grp_id = repmat(NaN,HDR.NS,1); 

	while ~isempty(tline),
		[n,v,sa]=str2double(tline);
		K = K+1; 
		HDR.FLT.channels.num(K,:)=n; 
			
		ch = n(1)+1;
		HDR.FLT.channels.seq(ch) = n(1);
		HDR.FLT.channels.id(ch) = n(2);
		HDR.FLT.channels.usage(ch) = n(3);
		HDR.Label{ch} = sa{4}; 
		HDR.FLT.channels.cal(ch) = n(5);
		HDR.FLT.channels.grd_mode(ch) = n(6);
		HDR.FLT.channels.grd_mode_name{ch} = sa{7};
		HDR.FLT.channels.grp_id(ch) = n(8);
		HDR.FLT.channels.no_sensors(ch) = n(9);

		for k=1:n(9);
			[tline,tch] = strtok(tch,[10,13]); 
			[n1,v1,sa]=str2double(tline);
			%sen = n1(1);
			sen = find(n1(1)==HDR.FLT.sensors.id);
			HDR.FLT.channels.Cal(ch,sen) = n1(2);
			if (length(sa)>=4) && ~strcmp(sa{4},HDR.FLT.sensors.name{sen}),
				fprintf(HDR.FILE.stderr,'Warning SOPEN(ET-MEG): sensor name does not fit: %s %s. \n    Maybe header of file %s is corrupted!\n',sa{4},HDR.FLT.sensors.name{sen}, HDR.FileName);
			end;
		end; 	
		[tline,tch] = strtok(tch,[10,13]);
	end; 	

	if isfield(HDR.FLT.Header,'name_of_data_file');
	if exist(fullfile(HDR.FILE.Path,HDR.FLT.Header.name_of_data_file)),
		HDR.FLT.datafile = HDR.FLT.Header.name_of_data_file; 
	end;
	end; 
	if HDR.FLT.Dataformat.type(1)<10; 
		HDR.Endianity = 'ieee-be'; 
	else	
		HDR.Endianity = 'ieee-le'; 
	end; 	
	switch mod(HDR.FLT.Dataformat.type(1),10)
	case 1,
		HDR.GDFTYP = 1;
	case 2,
		HDR.GDFTYP = 3;
	case 3,
		HDR.GDFTYP = 5;
	case 4,
		HDR.GDFTYP = 16;
	case 5,
		HDR.GDFTYP = 17;
	otherwise
		fprintf(HDR.FILE.stderr,'Error SOPEN(FLT): type %i not supported',type);
	end; 	

	[tmp,scale] = physicalunits(HDR.FLT.Measurement.sampling_unit);
	HDR.SampleRate = 10.^-HDR.FLT.Measurement.sampling_exponent ./ (HDR.FLT.Measurement.sampling_step.*scale);
	tmp = [HDR.FLT.Measurement.measurement_day,' ',HDR.FLT.Measurement.measurement_time];
	tmp (tmp=='.')=' ';
	tmp (tmp==':')=' ';
	HDR.T0([3,2,1,4:6]) = str2double(tmp);
	if isfield(HDR.FLT,'Patient'),
        	HDR.Patient.Sex = 0; % unknown
		HDR.FLT.Patient.remark = 'Do Not Modify !!!'; 
		if isfield(HDR.FLT.Patient,'sex')
        		tmp = deblank(HDR.FLT.Patient.sex);
	        	if length(tmp)
	        		HDR.Patient.Sex	= any(tmp(1)=='mM1') + 2*any(tmp(1)=='fFwW2');
	        	end;
	        end; 	
		if isfield(HDR.FLT.Patient,'birthday')
        		tmp = deblank(HDR.FLT.Patient.birthday);
        		tmp(tmp=='.' | tmp=='-' | tmp=='/')=' ';
        		[tmp,v,sa] = str2double(tmp); 
        		if length(tmp)==3 & ~any(v),
        			HDR.Patient.Birthday(1:3) = tmp;
        		end
        	end;
	end
	HDR.PhysDim = PhysDim_Group(HDR.FLT.channels.grp_id+1); 
	HDR.Cal = Cal_Group(HDR.FLT.channels.grp_id+1); 

	ix1 = strmatch('gr1s',HDR.FLT.channels.grd_mode_name,'exact');
	ix2 = find(any(HDR.FLT.channels.Cal(ix1,:),1));
	tmp = size(HDR.FLT.channels.Cal);
	ix1 = 1:tmp(1); ix2 = 1:tmp(2);
	grad.tra = HDR.FLT.channels.Cal;
	grad.pnt = HDR.FLT.sensors.XYZabcArea(:,1:3);
	grad.ori = HDR.FLT.sensors.XYZabcArea(:,4:6);
	grad.label = HDR.Label;
	grad.group_id = HDR.FLT.channels.grp_id;
	HDR.MEG.grad = grad; 

	%% compute HDR.ELEC.XYZ
	cal = HDR.FLT.channels.Cal;
	HDR.ELEC.XYZ = abs(cal)*(HDR.FLT.sensors.XYZabcArea(:,1:3));
	 
	HDR.Calib = sparse(2:HDR.NS+1,1:HDR.NS,HDR.Cal);
	if HDR.GDFTYP < 10,
		FLAG = HDR.FLAG; % backup
		HDR.FLAG.UCAL = 1; % default: no calibration information 
		% read scaling information
		ix  = strfind([HDR.FILE.Name,'.'],'.');
		fid = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.calib.txt']),'r'); 
		if fid < 0,
			fid = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name(1:ix(1)-1),'.calib.txt']),'r'); 
		end; 
		if fid>0,
			c   = fread(fid,[1,inf],'uint8=>char'); fclose(fid);
			[n1,v1,sa1] = str2double(c,9,[10,13]);
			%fid = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name(1:ix-1),'.calib2.txt']),'r');
			%c  = fread(fid,[1,inf],'uint8'); fclose(fid); 
			%[n2,v2,sa2] = str2double(c); 
			%HDR.Calib = sparse(2:HDR.NS+1,1:HDR.NS,(n2(4:131,5).*n1(4:131,3)./n2(4:131,4))./100);
			HDR.Calib = sparse(2:HDR.NS+1,1:HDR.NS,n1(4:131,3)*1e-13*(2^-12));
			HDR.FLAG = FLAG; % restore
		end
	end
	if ~isfield(HDR,'PhysDim') & ~isfield(HDR,'PhysDimCode') 
		HDR.PhysDimCode = zeros(HDR.NS,1); 
	end; 	
		
	% read data
	HDR.FILE.FID = -1; 
	if isfield(HDR.FLT,'datafile')
		fn = fullfile(HDR.FILE.Path,HDR.FLT.datafile);
		HDR.FILE.FID = fopen(fn,'rb',HDR.Endianity);
		if HDR.FILE.FID<0
			[p,f,e]=fileparts(HDR.FLT.datafile);
			fn = fullfile(HDR.FILE.Path,f);
			HDR.FILE.FID = fopen(fn,'rb',HDR.Endianity);
		end; 	
	else
		fn = fullfile(HDR.FILE.Path,HDR.FILE.Name);
		HDR.FILE.FID = fopen(fn,'rb',HDR.Endianity);
	end; 	

	if HDR.FILE.FID>0,
		HDR.FILE.OPEN=1; 
		fseek(HDR.FILE.FID,0,'eof'); 
		HDR.FILE.size = ftell(HDR.FILE.FID); 
		fseek(HDR.FILE.FID,0,'bof'); 
	else 
		fprintf(HDR.FILE.stderr,'Warning SOPEN(FLT): Binary data file %s not found!\n',fn); 		
	end;
	HDR.HeadLen  = 0; 
	HDR.FILE.POS = 0; 

        [datatyp,limits,datatypes,numbits,GDFTYP]=gdfdatatype(HDR.GDFTYP);
	HDR.AS.bpb    = HDR.NS*numbits/8;
	HDR.AS.endpos = HDR.FILE.size/HDR.AS.bpb;
        
	% check file size
        if (HDR.AS.bpb*HDR.NRec*HDR.SPR) ~= HDR.FILE.size,
        	% heuristic: try different data type
        	HDR.GDFTYP = 3; 
	        [datatyp,limits,datatypes,numbits,GDFTYP]=gdfdatatype(HDR.GDFTYP);
		HDR.AS.bpb    = HDR.NS*numbits/8;
		HDR.AS.endpos = HDR.FILE.size/HDR.AS.bpb;
        end;
        if (HDR.AS.bpb*HDR.NRec*HDR.SPR) ~= HDR.FILE.size,
        	% heuristic: try different data type
        	HDR.GDFTYP = 16; 
	        [datatyp,limits,datatypes,numbits,GDFTYP]=gdfdatatype(HDR.GDFTYP);
		HDR.AS.bpb    = HDR.NS*numbits/8;
		HDR.AS.endpos = HDR.FILE.size/HDR.AS.bpb;
        end;
        if (HDR.AS.bpb*HDR.NRec*HDR.SPR) ~= HDR.FILE.size,
        	fprintf(HDR.FILE.stderr,'Warning SOPEN(ET-MEG): size of file does not fit to header information\n');
        	fprintf(HDR.FILE.stderr,'\tFile:\t%s\n',fullfile(HDR.FILE.Path,HDR.FILE.Name));
        	fprintf(HDR.FILE.stderr,'\tFilesize:\t%i is not %i bytes\n',HDR.FILE.size,HDR.AS.bpb*HDR.NRec*HDR.SPR);
        	fprintf(HDR.FILE.stderr,'\tSamples:\t%i\n',HDR.NRec*HDR.SPR);
        	fprintf(HDR.FILE.stderr,'\tChannels:\t%i\n',HDR.NS);
        	fprintf(HDR.FILE.stderr,'\tDatatype:\t%s\n',datatyp);
        	HDR.SPR = floor(HDR.AS.endpos/HDR.NRec);
        end;	

elseif any(HDR.FILE.PERMISSION=='w'),
	fid = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.hdr']),'wt');
	if 0, isfield(HDR,'H1') 
		% copy header
		fwrite(fid,HDR.H1,'uint8');
	else	

		%%%%%%%% CHANNEL DATA %%%%%%%%%%%%%%%
		if ~isfield(HDR,'AS')
			HDR.AS.SampleRate = repmat(HDR.SampleRate,HDR.NS,1); 
		end;
		if ~isfield(HDR.AS,'SampleRate')
			HDR.AS.SampleRate = repmat(HDR.SampleRate,HDR.NS,1); 
		end;
		if ~isfield(HDR,'THRESHOLD')
			HDR.THRESHOLD = repmat(NaN,HDR.NS,2); 
		end;
		if ~isfield(HDR.Filter,'Notch')
			HDR.Filter.Notch = repmat(NaN,HDR.NS,1); 
		end;
		if ~isfield(HDR,'PhysDimCode')
			HDR.PhysDimCode = physicalunits(HDR.PhysDim); 
		end;
		if ~isfield(HDR,'LeadIdCode')
			HDR = leadidcodexyz(HDR); 
		end;
		if ~isfield(HDR,'REC')
			HDR.REC.Impedance = repmat(NaN,HDR.NS,1); 
		end;
		if ~isfield(HDR.REC,'Impedance')
			HDR.REC.Impedance = repmat(NaN,HDR.NS,1); 
		end;
		if ~isfield(HDR,'Off')
			HDR.Off = zeros(HDR.NS,1); 
		end;
		if ~isfield(HDR,'Cal')
			HDR.Cal = ones(HDR.NS,1); 
			HDR.Cal(HDR.InChanSelect) = diag(HDR.Calib(2:end,:));
		end;
		if length(HDR.Filter.HighPass)==1,
			HDR.Filter.HighPass = repmat(HDR.Filter.HighPass,HDR.NS,1); 
		end;
		if length(HDR.Cal)==1,
			HDR.Cal = repmat(HDR.Cal,HDR.NS,1); 
		end;
		if length(HDR.Filter.LowPass)==1,
			HDR.Filter.LowPass = repmat(HDR.Filter.LowPass,HDR.NS,1); 
		end;
		if length(HDR.Filter.Notch)==1,
			HDR.Filter.Notch = repmat(HDR.Filter.Notch,HDR.NS,1); 
		end;

		PhysDim = physicalunits(HDR.PhysDimCode); 
		if ~isfield(HDR,'FLT'); HDR.FLT = []; end; 

		%%%%%%%%% FIXED HEADER %%%%%%%%%%%%%%		
		fprintf(fid,'[Header]\n'); 
		fprintf(fid,'version=2.1\nid=1\n'); 
		fprintf(fid,'name=created by BioSig for Octave and Matlab http://biosig.sf.net/\n'); 
		fprintf(fid,'comment= -\n'); 
		fprintf(fid,'name_of_data_file=%s\n',HDR.FILE.Name); 

		fprintf(fid,'\n[Dataformat]\n'); 
		fprintf(fid,'* data types : HP-UX data\n');
		fprintf(fid,'*               1=(1 byte int)    2=(2 byte int)    3=(4 byte int)\n');
		fprintf(fid,'*               4=(4 byte float)  5=(8 byte float)  6=(ASCII)\n');
		fprintf(fid,'*              LINUX data\n');
		fprintf(fid,'*              11=(1 byte int)   12=(2 byte int)   13=(4 byte int)\n');
		fprintf(fid,'*              14=(4 byte float) 15=(8 byte float) 16=(ASCII)\n');
		fprintf(fid,'version=1.0\nid=1\nname=ET-MEG double data format\n'); 
		fprintf(fid,'type=14\n'); HDR.GDFTYP=16; % float32
		if isfield(HDR,'data')
			[HDR.SPR,HDR.NS]=size(HDR.data);
			HDR.NRec = 1;
		end;	
		fprintf(fid,'number_of_samples=%i\n',HDR.NRec*HDR.SPR); 

		fprintf(fid,'\n[Measurement]\nversion=0.0\nlaboratory_name=%s\n',HDR.FLT.Measurement.laboratory_name); 
		fprintf(fid,'*\nmeasurement_day=%02i.%02i.%04i\n',HDR.T0([3,2,1])); 
		fprintf(fid,'measurement_time=%02i:%02i:%02i\n',HDR.T0(4:6));
		fprintf(fid,'sampling_unit=s\n'); 
		e = floor(log10(1/HDR.SampleRate));
		fprintf(fid,'sampling_exponent=%i\n',e); 
		fprintf(fid,'sampling_step=%f\n',10^-e/HDR.SampleRate); 
		 
		if isfield(HDR,'Patient')
			fprintf(fid,'\n[Patient]\nid= 1\nversion= 1.00\n'); 
			fprintf(fid,'name= -\n');   %%% will not be supported due to privacy reasons
			if HDR.Patient.Sex==1,
				fprintf(fid,'sex= male\n');   
			elseif HDR.Patient.Sex==2,
				fprintf(fid,'sex= female\n');   
			else
				fprintf(fid,'sex= -\n');
			end;	
			if isfield(HDR.Patient,'Birthday')
				fprintf(fid,'birthday=%02i.%02i.%04i\n',HDR.T0([3,2,1]));
			else 	
				fprintf(fid,'birthday= -\n');
			end; 	
		end; 	 
		 
		fprintf(fid,'\n[System]\nversion=0.0\n'); 
		fprintf(fid,'number_of_channels=%i\n',HDR.NS); 
		%fprintf(fid,'SamplingRate=%i\n',HDR.SampleRate); 


		fprintf(fid,'*---------------------------------------------------------------\n');
		fprintf(fid,'*seq id   u name              calib grd grd_name  grp  n_sensors\n');
		fprintf(fid,'*---------------------------------------------------------------\n');
		fprintf(fid,'parameter_of_channels={\n');
		if isfield(HDR.FLT,'channels')
			% write original channel header 
			num = [HDR.FLT.channels.seq,HDR.FLT.channels.id,HDR.FLT.channels.usage,HDR.FLT.channels.cal,HDR.FLT.channels.grd_mode,HDR.FLT.channels.grp_id];       
			for k=1:HDR.NS,
				ix = find(HDR.FLT.channels.Cal(k,:));
				fprintf(fid,'%04i %04i %i %-17s %5.3f  %i  %-8s  %04i  %i\n',num(k,1:3),HDR.Label{k},num(k,4:5),HDR.FLT.channels.grd_mode_name{k},num(k,6),length(ix)); 
				for k1 = 1:length(ix)
					fprintf(fid,'\t%04i %9.6f * %s\n',ix(k1)-1,full(HDR.FLT.channels.Cal(k,ix(k1))),HDR.FLT.sensors.name{ix(k1)});
				end; 
			end; 	
		else
			for k=1:HDR.NS,
				fprintf(fid,'%04i  %04i 1 %-16s 1.000  1  ####      0000   1\n',k-1,k-1,HDR.Label{k});
				fprintf(fid,'        %04i 1.000000 * CH%i\n',k-1,k-1);
			end; 	
		end;
		
		if isfield(HDR.FLT,'sensors'); 
			fprintf(fid,'}\n\nnumber_of_sensors=%i\n',HDR.FLT.System.number_of_sensors);
		else
			fprintf(fid,'}\n\nnumber_of_sensors=%i\n',HDR.NS);
		end;	
		fprintf(fid,'*----------------------------------------------------------------------------------------\n');
		fprintf(fid,'*id  name     type mod    x         y         z         a         b         c        area\n');
		fprintf(fid,'*----------------------------------------------------------------------------------------\n');
		fprintf(fid,'parameter_of_sensors={\n');
		if isfield(HDR.FLT,'sensors'); 
			for k=1:size(HDR.FLT.sensors.id,1),
				fprintf(fid,'%04i %-10s %i %04i  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n',HDR.FLT.sensors.id(k),HDR.FLT.sensors.name{k},HDR.FLT.sensors.type(k),HDR.FLT.sensors.mod(k),HDR.FLT.sensors.XYZabcArea(k,:)); 
			end; 	
		else
			for k=1:HDR.NS,
				fprintf(fid,'%04i CH%-08i  1 0000  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %11.9f\n',k-1,k-1,HDR.ELEC.XYZ(k,1:3),[NaN,NaN,NaN,NaN]);
			end;
		end;
		 	
		fprintf(fid,'}\n\nnumber_of_groups=%i\n',HDR.FLT.System.number_of_groups);
		fprintf(fid,'*----------------------------------------\n');
		fprintf(fid,'*id  u name             unit   exp  calib\n');
		fprintf(fid,'*----------------------------------------\n');
		fprintf(fid,'parameter_of_groups={\n');
		if isfield(HDR.FLT,'groups')
			for k=1:size(HDR.FLT.groups.id,1),
				fprintf(fid,'%04i %i %-16s %-6s %i %8.3f\n',HDR.FLT.groups.id(k),HDR.FLT.groups.usage(k),HDR.FLT.groups.name{k},HDR.FLT.groups.unit{k},HDR.FLT.groups.unit_exp(k),HDR.FLT.groups.calib(k)); 
			end; 	
		else
			fprintf(fid,'0001 1 ET-Mag_80WH      T      0    1.000\n');
			fprintf(fid,'0002 1 ET-AxGrd_80WH    T      0    1.000\n');
			fprintf(fid,'0003 1 ET-PlGrd_80WH    T      0    1.000\n');
			fprintf(fid,'0004 1 ET-Mag_RefCh     T      0    1.000\n');
			fprintf(fid,'0005 1 ET-AxGrd_RefCh   T      0    1.000\n');
			fprintf(fid,'0006 1 ET-PlGrd_RefCh   T      0    1.000\n');
			fprintf(fid,'0007 1 Trigger          V      0    1.000\n');
			fprintf(fid,'0008 1 EEG              V      0    1.000\n');
			fprintf(fid,'0009 1 ECG              V      0    1.000\n');
			fprintf(fid,'0010 1 Etc              V      0    1.000\n');
			fprintf(fid,'0011 0 Null_Channel     V      0    1.000\n');
		end;	
		fprintf(fid,'}\n');
		fprintf(fid,'\nnumber_of_modules=%i\n',HDR.FLT.System.number_of_modules);
		fprintf(fid,'*-------------------------------------------------------------------------\n');
		fprintf(fid,'*id  name          x      y      z      a      b      c      unit exp name\n');
		fprintf(fid,'*-------------------------------------------------------------------------\n');
		fprintf(fid,'parameter_of_modules={\n');
		if isfield(HDR.FLT,'modules')
			for k=1:size(HDR.FLT.modules.id,1),
				fprintf(fid,'%04i %-13s %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f %i %s\n',HDR.FLT.modules.id(k),HDR.FLT.modules.name{k},HDR.FLT.modules.XYZabc(k,:),HDR.FLT.modules.unit(k),HDR.FLT.modules.exp(k),HDR.FLT.modules.unitname{k}); 
			end; 	
			fprintf(fid,'}\n'); 
		else			fprintf(fid,'0000 Magnetic      0.000  0.000  0.000  0.000  0.000  0.000  1.000 0 m\n');
			fprintf(fid,'0007 Electric      0.000  0.000  0.000  0.000  0.000  0.000  1.000 0 m\n}\n');
		end;
	end;
	fclose(fid);

	HDR.FILE.FID  = fopen(fullfile(HDR.FILE.Path,HDR.FILE.Name),'wb','ieee-le');
	HDR.FILE.OPEN = 2;
	HDR.HeadLen   = 0; 


elseif any(HDR.FILE.PERMISSION=='a'),
	HDR.FILE.FID  = fopen(fullfile(HDR.FILE.Path,HDR.FILE.Name),'ab','ieee-le');
	HDR.FILE.OPEN = 3;
	HDR.HeadLen   = 0; 

end;

