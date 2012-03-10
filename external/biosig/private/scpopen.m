function [HDR]=scpopen(arg1,CHAN,arg4,arg5,arg6)
% SCPOPEN reads and writes SCP-ECG files 
%
% SCPOPEN is an auxillary function to SOPEN for 
% opening of SCP-ECG files for reading ECG waveform data
% 
% Use SOPEN instead of SCPOPEN  
% 
% See also: fopen, SOPEN, 


%	$Id$
%	(C) 2004,2006,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BioSig-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

if nargin<2, CHAN=0; end;

if isstruct(arg1) 
        HDR=arg1; 
        FILENAME=HDR.FileName;
elseif ischar(arg1); 
        HDR.FileName=arg1;
        fprintf(2,'Warning SCPOPEN: the use of SCPOPEN is discouraged; please use SOPEN instead.\n');
end;

VER = version;

fid = fopen(HDR.FileName,HDR.FILE.PERMISSION,'ieee-le');
HDR.FILE.FID = fid; 
if ~isempty(findstr(HDR.FILE.PERMISSION,'r')),		%%%%% READ 
	tmpbytes = fread(fid,inf,'uchar');
        tmpcrc   = crc16eval(tmpbytes(3:end));
	fseek(fid, 0, 'bof'); 
        HDR.FILE.CRC = fread(fid,1,'uint16');
	if (HDR.FILE.CRC ~= tmpcrc);
		fprintf(HDR.FILE.stderr,'Warning: CRC check failed (%x vs %x)\n',tmpcrc,HDR.FILE.CRC);        
        end;
	
        HDR.FILE.Length = fread(fid,1,'uint32');
        if HDR.FILE.Length~=HDR.FILE.size,
                fprintf(HDR.FILE.stderr,'Warning SCPOPEN: header information contains incorrect file size %i %i \n',HDR.FILE.Length,HDR.FILE.size);
        end; 
	HDR.data = [];
        
        DHT = [0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,7,-7,8,-8,9,-9;0,1,5,3,11,7,23,15,47,31,95,63,191,127,383,255,767,511,1023]';
        prefix  = [1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,10,10];
	PrefixLength = prefix; 

        %PREFIX = [0,4,5,12,13,28,29,60,61,124,125,252,253,508,509,1020,1021,1022,1023];
        PREFIX  = [0,4,5,12,13,28,29,60,61,124,125,252,253,508,509,1020,1021,1022,1023]'.*2.^[32-prefix]';
        codelength = [1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,18,26];
        mask    = [1,7,7,15,15,31,31,63,63,127,127,255,255,511,511,1023,1023,1023,1023]'.*2.^[32-prefix]';
        %MASK    = dec2bin(mask);
        %mask   = [1,7,7,15,15,31,31,63,63,127,127,255,255,511,511,1023,1023,1023,1023]';

        mask2    = [1,7,7,15,15,31,31,63,63,127,127,255,255,511,511,1023,1023,1023,1023]';
        PREFIX2  = DHT(:,2);

	HT19999 = [prefix',codelength',ones(length(prefix),1),DHT];
	HT = [prefix',codelength',ones(length(prefix),1),DHT];
        
        dd = [0:255]';
        ACC = zeros(size(dd));
        c = 0;
        for k2 = 1:8,
                ACC = ACC + (dd>127).*(2^c);
                dd  = mod(dd*2, 256);
                c   = c + 1;
        end;
        
        section.CRC     = fread(fid,1,'uint16');
        section.ID      = fread(fid,1,'uint16');
        section.Length  = fread(fid,1,'uint32');
        section.Version = fread(fid,[1,2],'uint8');
        section.tmp     = fread(fid,[1,6],'uint8');
        
        NSections = min(11,(section.Length-16)/10);
        for k = 1:NSections,
                HDR.Block(k).id = k; 
                HDR.Block(k).length = 0; 
                HDR.Block(k).startpos = -1;
        end;
	for K = 1:NSections,
                k = fread(fid,1,'uint16');
		len = fread(fid,1,'uint32');
		pos = fread(fid,1,'uint32');
		if ((k > 0) && (k < NSections))
	                HDR.Block(k).id = k; 
        	        HDR.Block(k).length = len; 
                	HDR.Block(k).startpos = pos-1;

%% [HDR.Block(k).id ,length(tmpbytes), HDR.Block(k).length, HDR.Block(k).length+HDR.Block(k).startpos]                
		%% FIXME: instead of min(...,FileSize) a warning or error message should be reported 
	                tmpcrc = crc16eval(tmpbytes(HDR.Block(k).startpos+3:min(HDR.Block(k).startpos+HDR.Block(k).length,HDR.FILE.size)));

	                if (HDR.Block(k).length>0)
        	        if (tmpcrc~=(tmpbytes(HDR.Block(k).startpos+(1:2))'*[1;256]))
                		fprintf(HDR.FILE.stderr,'Warning SCPOPEN: faulty CRC %04x in section %i\n',tmpcrc,k-1);
	                end;
        	        end;
		end;
        end;
        
%%[[HDR.Block.id];[HDR.Block.length];[HDR.Block.startpos]]'
	
	% default values - in case Section 6 is missing  
	HDR.NS = 0; HDR.SPR = 0; HDR.NRec = 0; HDR.Calib = zeros(1,0); 
        secList = find([HDR.Block.length]);
        for K = secList(1:end),
                if fseek(fid,HDR.Block(K).startpos,'bof');
                        fprintf(HDR.FILE.stderr,'Warning SCPOPEN: section %i not available, although it is listed in Section 0\n',secList(K+1));
                end;
                section.CRC     = fread(fid,1,'uint16');
                section.ID      = fread(fid,1,'uint16');
                section.Length  = fread(fid,1,'uint32');
                section.Version = fread(fid,[1,2],'uint8');
                section.tmp     = fread(fid,[1,6],'uint8'); 
	
                HDR.SCP.Section{find(K==secList)} = section;
                if (section.Length==0),
                elseif section.ID==0, 
                        NSections = (section.Length-16)/10;
                        for k = 1:NSections,
                                HDR.Block(k).id = fread(fid,1,'uint16');    
                                HDR.Block(k).length = fread(fid,1,'uint32');    
                                HDR.Block(k).startpos = fread(fid,1,'uint32')-1;    
                        end;
                        
                elseif section.ID==1,
                        tag = 0; 
                        k1  = 0;
                        Sect1Len = section.Length-16;
                        ListOfRequiredTags = [2,14,25,26];
                        ListOfRecommendedTags = [0,1,5,8,15,34];
                        while (tag~=255) & (Sect1Len>2),
                                tag = fread(fid,1,'uint8');
                                len = fread(fid,1,'uint16');
                                Sect1Len = Sect1Len - 3 - len; 
%% [tag,len,Sect1Len],          %% DEBUGGING information
                                field = fread(fid,[1,len],'uchar');
                                if tag == 0,	
                                        ListOfRecommendedTags(ListOfRecommendedTags==tag)=[];
                                        HDR.Patient.Name = char(field);  %% LastName
                                elseif tag == 1,
                                        ListOfRecommendedTags(ListOfRecommendedTags==tag)=[];
                                        HDR.Patient.FirstName = char(field);
                                elseif tag == 2,
                                        ListOfRequiredTags(find(ListOfRequiredTags==2))=[];
                                        HDR.Patient.Id = char(field);
                                elseif tag == 3,
                                        HDR.Patient.LastName2 = char(field);
                                elseif tag == 4,
                                        HDR.Patient.Age = field(1:2)*[1;256];
                                        tmp = field(3);
                                        if     tmp==1, HDR.Patient.Age = HDR.Patient.Age; % unit='Y';
                                        elseif tmp==2, HDR.Patient.Age = HDR.Patient.Age/12; % unit='M';
                                        elseif tmp==3, HDR.Patient.Age = HDR.Patient.Age/52; % unit='W';
                                        elseif tmp==4, HDR.Patient.Age = HDR.Patient.Age/365.25; % unit='d';
                                        elseif tmp==5, HDR.Patient.Age = HDR.Patient.Age/(365.25*24); %unit='h';
                                        else warning('units of age not specified');
                                        end;
                                elseif (tag == 5) 
                                        ListOfRecommendedTags(ListOfRecommendedTags==tag)=[];
                                	if any(field(1:4)~=0)
                                        	HDR.Patient.Birthday = [field(1:2)*[1;256],field(3:4),12,0,0];
                                        end;	
                                elseif (tag == 6) 
                                	if any(field(1:3)),
                                        HDR.Patient.Height = field(1:2)*[1;256];
                                        tmp = field(3);
                                        if tmp==1, % unit='cm';
                                        elseif tmp==2, HDR.Patient.Height = HDR.Patient.Height*2.54; %unit='inches'; 
                                        elseif tmp==3, HDR.Patient.Height = HDR.Patient.Height*0.1; %unit='mm';
                                        else warning('units of height not specified');
                                        end;
                                        end;
                                elseif (tag == 7) 
                                	if any(field(1:3)),
                                        HDR.Patient.Weight = field(1:2)*[1;256];
                                        tmp = field(3);
                                        if tmp==1, % unit='kg';
                                        elseif tmp==2, HDR.Patient.Weight = HDR.Patient.Weight/1000; %unit='g';
                                        elseif tmp==3, HDR.Patient.Weight = HDR.Patient.Weight/2.2; %unit='pound';
                                        elseif tmp==4, HDR.Patient.Weight = HDR.Patient.Weight*0.0284; %unit='ounce';
                                        else warning('units of weight not specified');
                                        end;
                                        end;
                                elseif tag == 8,
                                        ListOfRecommendedTags(ListOfRecommendedTags==tag)=[];
                                        HDR.Patient.Sex = field;
                                elseif tag == 9,
                                        HDR.Patient.Race = field;
                                elseif tag == 10,
					if (field(1)~=0)
	                                        HDR.Patient.Medication = field;
					else	
	                                        HDR.Patient.Medication.Code = field(2:3);
						HDR.Patient.Medication = field(4:end);
					end;	
                                elseif tag == 11,
                                        HDR.Patient.BloodPressure.Systolic = field*[1;256];
                                elseif tag == 12,
                                        HDR.Patient.BloodPressure.Diastolic = field*[1;256];
                                elseif tag == 13,
                                        HDR.Patient.Diagnosis = char(field);
                                elseif tag == 14,
                                        ListOfRequiredTags(ListOfRequiredTags==tag)=[];
                                        HDR.SCP1.AcquiringDeviceID = char(field);
                                        HDR.VERSION = field(15)/10;
                                elseif tag == 15,
                                        ListOfRecommendedTags(ListOfRecommendedTags==tag)=[];
                                        HDR.SCP1.AnalyisingDeviceID = char(field);
                                elseif tag == 16,
                                        HDR.SCP1.AcquiringInstitution = char(field);
                                elseif tag == 17,
                                        HDR.SCP1.AnalyzingInstitution = char(field);
                                elseif tag == 18,
                                        HDR.SCP1.AcquiringDepartment = char(field);
                                elseif tag == 19,
                                        HDR.SCP1.AnalyisingDepartment = char(field);
                                elseif tag == 20,
                                        HDR.SCP1.Physician = char(field);
                                elseif tag == 21,
                                        HDR.SCP1.LatestComfirmingPhysician = char(field);
                                elseif tag == 22,
                                        HDR.SCP1.Technician = char(field);
                                elseif tag == 23,
                                        HDR.SCP1.Room = char(field);
                                elseif tag == 24,
                                        HDR.SCP1.Emergency = field;
                                elseif tag == 25,
                                        ListOfRequiredTags(ListOfRequiredTags==tag)=[];
                                        HDR.T0(1,1:3) = [field(1:2)*[1;256],field(3:4)];
                                elseif tag == 26,
                                        ListOfRequiredTags(ListOfRequiredTags==tag)=[];
                                        HDR.T0(1,4:6) = field(1:3);
                                elseif tag == 27,
                                        HDR.Filter.HighPass = field(1:2)*[1;256]/100;
                                elseif tag == 28,
                                        HDR.Filter.LowPass = field(1:2)*[1;256]/100;
                                elseif tag == 29,
                                        if (field==0)
                                                HDR.FILTER.Notch = NaN; 
                                        elseif bitand(field,1)
                                                HDR.FILTER.Notch = 60; % 60Hz Notch 
                                        elseif bitand(field,2)
                                                HDR.FILTER.Notch = 50; % 50Hz Notch 
                                        elseif bitand(field,3)==0
                                                HDR.FILTER.Notch = -1; % Notch Off
                                        end;
                                        HDR.SCP1.Filter.BitMap = field;
                                elseif tag == 30,
                                        HDR.SCP1.FreeText = char(field);
                                elseif tag == 31,
                                        HDR.SCP1.ECGSequenceNumber = char(field);
                                elseif tag == 32,
                                        HDR.SCP1.MedicalHistoryCodes = char(field);
                                elseif tag == 33,
                                        HDR.SCP1.ElectrodeConfigurationCodes = field;
                                elseif tag == 34,
                                        ListOfRecommendedTags(ListOfRecommendedTags==tag)=[];
                                        HDR.SCP1.Timezone = field;
                                elseif tag == 35,
                                        HDR.SCP1.MedicalHistory = char(field);
                                elseif tag == 255,
                                        % section terminator	
                                elseif tag >= 200,
                                	% manufacturer specific - not standardized 
                                else
                                        fprintf(HDR.FILE.stderr,'Warning SCOPEN: unknown tag %i (section 1)\n',tag);
                                end;
                        end;
                        if ~isempty(ListOfRequiredTags)
                                fprintf(HDR.FILE.stderr,'Warning SCPOPEN: the following tags are required but missing in file %s\n',HDR.FileName);
                                disp(ListOfRequiredTags);
                        end;
                        if ~isempty(ListOfRecommendedTags)
                                fprintf(HDR.FILE.stderr,'Warning SCPOPEN: the following tags are recommended but missing in file %s\n',HDR.FileName);
                                disp(ListOfRecommendedTags);
                        end;
                        
                elseif section.ID==2, 	% Huffman tables 
                        HDR.SCP2.NHT = fread(fid,1,'uint16');            
                        HDR.SCP2.NCT = fread(fid,1,'uint16');    
			if HDR.SCP2.NHT~=19999,
				NHT = HDR.SCP2.NHT;
			else
				NHT = 0; 
			end;
                        k3 = 0;
                        for k1 = 1:NHT,
                        	HT1 = zeros(HDR.SCP2.NCT,5);
                                for k2 = 1:HDR.SCP2.NCT,
                                	tmp = fread(fid,3,'uint8') ;
                                        HDR.SCP2.prefix = tmp(1);	% PrefixLength
                                        HDR.SCP2.codelength = tmp(2);	% CodeLength
                                        HDR.SCP2.TableModeSwitch = tmp(3);	
                                        tmp(4) = fread(fid,1,'int16');  % BaseValue   
                                        tmp(5) = fread(fid,1,'uint32'); % BaseCode    
                                	k3 = k3   + 1;
				        HT (k3,:) = [tmp']; 
				        HT1(k2,:) = [tmp']; 
                                end;
                                HDR.SCP2.HTree{k1} = makeTree(HT1);
                                HDR.SCP2.HTs{k1} = HT1;
			end;
			if HDR.SCP2.NHT~=19999,
				HDR.SCP2.HT = HT;
			else
				tmp = size(HT19999,1);
				HDR.SCP2.HT = [ones(tmp,1),[1:tmp]',HT19999];
                                HDR.SCP2.HTree{1} = makeTree(HT19999);
                                HDR.SCP2.HTs{1} = HT19999;
			end;

                elseif section.ID==3, 
                        HDR.NS = fread(fid,1,'uint8');
                        HDR.FLAG.Byte = fread(fid,1,'uint8');    
                        if ~bitand(HDR.FLAG.Byte,4)
                                fprintf(HDR.FILE.stdout,'Warning SCPOPEN: not all leads simultaneously recorded - this mode is not supported.\n');
                        end;
                                        
                        HDR.FLAG.ReferenceBeat = mod(HDR.FLAG.Byte,2);    
                        %HDR.NS = floor(mod(HDR.FLAG.Byte,128)/8);    
                        for k = 1:HDR.NS,
                                HDR.LeadPos(k,1:2) = fread(fid,[1,2],'uint32');    
                                HDR.LeadIdCode(k,1) = fread(fid,1,'uint8');    
                        end;
                        HDR.N = max(HDR.LeadPos(:))-min(HDR.LeadPos(:))+1;
                        HDR.AS.SPR = HDR.LeadPos(:,2)-HDR.LeadPos(:,1)+1;
                        HDR.SPR = HDR.AS.SPR(1);  
			for k = 2:HDR.NS
				HDR.SPR = lcm(HDR.SPR,HDR.AS.SPR(k)); 
			end; 	                        
                        
                        HDR = leadidcodexyz(HDR);
                        for k = 1:HDR.NS,
                                if 0,
                                elseif (HDR.LeadIdCode(k)==0),
                                        HDR.Label{k} = 'unspecified lead';
                                elseif (HDR.VERSION <= 1.3) & (HDR.LeadIdCode(k) < 86),
                                %        HDR.Label{k} = H.Label(H.LeadIdCode==HDR.LeadIdCode(k));
                                elseif (HDR.VERSION <= 1.3) & (HDR.LeadIdCode(k) > 99),
                                        HDR.Label{k} = 'manufacturer specific';
                                elseif (HDR.VERSION >= 2.0) & (HDR.LeadIdCode(k) < 151),
                                %        HDR.Label{k} = H.Label(H.LeadIdCode==HDR.LeadIdCode(k));
                                elseif (HDR.VERSION >= 2.0) & (HDR.LeadIdCode(k) > 199),
                                        HDR.Label{k} = 'manufacturer specific';
                                else
                                        HDR.Label{k} = 'reserved';
                                end;
                        end;
                        HDR.Label = strvcat(HDR.Label);

                elseif section.ID==4, 
                        HDR.SCP4.L = fread(fid,1,'int16');    
                        HDR.SCP4.fc0 = fread(fid,1,'int16');    
                        HDR.SCP4.N = fread(fid,1,'int16');    
                        HDR.SCP4.type = fread(fid,[7,HDR.SCP4.N],'uint16')'*[1,0,0,0; 0,1,0,0;0,2^16,0,0; 0,0,1,0;0,0,2^16,0; 0,0,0,1;0,0,0,2^16];   

                        tmp = fread(fid,[2*HDR.SCP4.N],'uint32');
                        HDR.SCP4.PA = reshape(tmp,2,HDR.SCP4.N)';   
                        HDR.SCP4.pa = [0;tmp;HDR.N];   
                        
                elseif any(section.ID==[5,6]), 

                        SCP = [];
                        SCP.Cal = fread(fid,1,'int16')/1e6;    % quant in nV, converted into mV
                        SCP.PhysDim = 'mV';
                        SCP.Dur = fread(fid,1,'int16');    
                        SCP.SampleRate = 1e6/SCP.Dur;
                        SCP.FLAG.DIFF  = fread(fid,1,'uint8');    
                        SCP.FLAG.bimodal_compression = fread(fid,1,'uint8');    

                        if isnan(HDR.NS),
				HDR.ERROR.status = -1; 
				HDR.ERROR.message = sprintf('Error SCPOPEN: could not read %s\n',HDR.FileName);
				fprintf(HDR.FILE.stderr,'Error SCPOPEN: could not read %s\n',HDR.FileName);
				return;
			end;
			
			if CHAN==0, CHAN = 1:HDR.NS; end;
                        SCP.SPR = fread(fid,HDR.NS,'uint16');
			HDR.InChanSelect = CHAN; 

                        if section.ID==6,
                                HDR.HeadLen = ftell(fid);
                                HDR.FLAG.DIFF = SCP.FLAG.DIFF;
                                HDR.FLAG.bimodal_compression = SCP.FLAG.bimodal_compression;
                                HDR.data = [];
                                outlen = HDR.SPR; 
			        HDR.Calib = sparse(2:HDR.NS+1, 1:HDR.NS, SCP.Cal);
                        elseif isfield(HDR,'SCP4') %% HACK: do no know whether it is correct  
                        	outlen = floor(1000*HDR.SCP4.L/SCP.Dur);
                        else 
                        	outlen = inf;
                        end;

                        if ~isfield(HDR,'SCP2'),
                                if any(SCP.SPR(1)~=SCP.SPR),
                                        error('SCPOPEN: SPR do not fit');
                                else
                                        S2 = fread(fid,[SCP.SPR(1)/2,HDR.NS],'int16');
                                end;
				%S2 = S2(:,CHAN); 
        
                        elseif (HDR.SCP2.NHT==1) && (HDR.SCP2.NCT==1) && (HDR.SCP2.prefix==0), 
				codelength = HDR.SCP2.HT(1,4);
                                if (codelength==16)
                                        S2 = fread(fid,[HDR.N,HDR.NS],'int16');  
                                elseif (codelength==8)
                                        S2 = fread(fid,[HDR.N,HDR.NS],'int8');  
                                else
                                        fprintf(HDR.FILE.stderr,'Warning SCPOPEN: codelength %i is not supported yet.',codelength);
                                        fprintf(HDR.FILE.stderr,' Contact <a.schloegl@ieee.org>\n');
                                        return;
                                end;
				%S2 = S2(:,CHAN); 
                                
                        elseif 1, HDR.SCP2.NHT~=19999;
                        	%% User specific Huffman table 
                        	%% a more elegant Huffman decoder is used here %%
                                for k = 1:HDR.NS,
                                        SCP.data{k} = fread(fid,SCP.SPR(k),'uint8');    
                                end;
%                                S2 = repmat(NaN,outlen,length(HDR.InChanSelect));
                                clear S2;  
                                sz = inf;
                                for k3 = 1:length(HDR.InChanSelect), k = HDR.InChanSelect(k3); %HDR.NS,
					outdata{k3} = DecodeHuffman(HDR.SCP2.HTree,HDR.SCP2.HTs,SCP.data{k},outlen);
					sz = min(sz,length(outdata{k3}));
				end;
				
                                for k3 = 1:length(HDR.InChanSelect), k = HDR.InChanSelect(k3); %HDR.NS,
                                	S2(:,k) = outdata{k3}(1:sz);
                                end; 	 
				accu=0;  	

                        elseif HDR.SCP2.NHT==19999,
                                HuffTab = DHT;
                                for k = 1:HDR.NS,
                                        SCP.data{k} = fread(fid,SCP.SPR(k),'uint8');    
                                end;
                                %for k = 1:HDR.NS,
                                for k3 = 1:length(HDR.InChanSelect), k = HDR.InChanSelect(k3); %HDR.NS,
                                %for k = CHAN(:)',
                                        s2 = SCP.data{k};
                                        s2 = [s2; repmat(0,ceil(max(HDR.SCP2.HT(:,4))/8),1)];
					k1 = 0;	
					l2 = 0; 
					accu = 0;
					c  = 0; 
					x  = [];
					HT = HDR.SCP2.HT(find(HDR.SCP2.HT(:,1)==1),3:7);
					while (l2 < HDR.LeadPos(k,2)),
						while ((c < max(HT(:,2))) & (k1<length(s2)-1));
							k1 = k1 + 1;
							dd = s2(k1);
							accu = accu + ACC(dd+1)*(2^c);
							c = c + 8;

							if 0, %for k2 = 1:8,
								accu = accu + (dd>127)*(2^c);
								dd = mod(dd*2,256);
								c = c + 1;
							end;
						end;

                                                ixx = 1;
                                                %acc = mod(accu,2^32);   % bitand returns NaN if accu >= 2^32
						acc = accu - 2^32*fix(accu*(2^(-32)));   % bitand returns NaN if accu >= 2^32
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
							%acc2  = fix(accu*(2^(-HT(ixx,1))));
							%tmp = mod(fix(accu*(2^(-HT(ixx,1)))),2^dd);
							
                                                        tmp = fix(accu*(2^(-HT(ixx,1))));       % bitshift(accu,-HT(ixx,1))
                                                        tmp = tmp - (2^dd)*fix(tmp*(2^(-dd)));  % bitand(...,2^dd)
                                                        
                                                        %tmp = bitand(accu,(2^dd-1)*(2^HT(ixx,1)))*(2^-HT(ixx,1));
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
						accu = fix(accu*2^(-HT(ixx,2)));
						c = c - HT(ixx,2); 
					end;
					x = x(:);
                                        if k3==1,
                                                S2=x(:,ones(1,k));
                                        elseif size(x,1)==size(S2,1),
                                                S2(:,k) = x;
					else
	                                        fprintf(HDR.FILE.stderr,'Error SCPOPEN: Huffman decoding failed (%i) \n',size(x,1));
	    					HDR.data = S2;
						return;
                                        end;
				end;
                                
                                
                        elseif (HDR.SCP2.NHT==19999), % alternative decoding algorithm. 
                                warning('this branch is experimental - it might be broken')
                                HuffTab = DHT;
                                for k = 1:HDR.NS,
                                        SCP.data{k} = fread(fid,SCP.SPR(k),'uint8');
                                end;
                                %for k = 1:HDR.NS,
                                for k3 = 1:length(HDR.InChanSelect), k = HDR.InChanSelect(k3); %HDR.NS,
                                %for k = CHAN(:)',
				        tmp = SCP.data{k};
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
                                                while (bitand(accu,mask(ixx)) ~= PREFIX(ixx)), 
                                                        ixx = ixx + 1;
                                                end;

                                                if ixx < 18,
                                                        c = c + prefix(ixx);
                                                        %accu  = bitshift(accu, prefix(ixx),32);
                                                        accu  = mod(accu.*(2^prefix(ixx)),2^32);
                                                        l2    = l2 + 1;
                                                        x(l2) = HuffTab(ixx,1);
                                                        
                                                elseif ixx == 18,
                                                        c = c + prefix(ixx) + 8;
                                                        %accu = bitshift(accu, prefix(ixx),32);
                                                        accu  = mod(accu.*(2^prefix(ixx)),2^32);
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
                                                        c = c + prefix(ixx);
                                                        %accu = bitshift(accu, prefix(ixx),32);
                                                        accu  = mod(accu.*(2^prefix(ixx)),2^32);
                                                        l2    = l2 + 1;
                                                        while (c > 7) & (l < Ntmp),
                                                                l = l+1;
                                                                c = c-8;
                                                                accu = accu + tmp(l)*2^c;
                                                        end;
                                                        
                                                        acc1 = mod(floor(accu*2^(-16)),2^16);
                                                        %accu = bitshift(accu, 16, 32);
                                                        accu = mod(accu.*(2^16), 2^32);
                                                        
                                                        x(l2) = acc1-(acc1>=2^15)*2^16;
                                                        acc2 = 0;
                                                        for kk=1:16,
                                                                acc2 = acc2*2+mod(acc1,2);
                                                                acc1 = floor(acc1/2);
                                                        end;
                                                        %x(l2) = acc2;
                                                        c = c + 16;
                                                end;
                                                
                                                while (c > 7) & (l < Ntmp),
                                                        l = l+1;
                                                        c = c-8;
                                                        accu = accu + tmp(l)*(2^c);
                                                end;
                                        end;

                                        x = x(1:end-1)';
                                        if k3==1,
                                                S2=x(:,ones(1,k));
                                        elseif size(x,1)==size(S2,1),
                                                S2(:,k) = x;
                                        elseif 1,
	                                        fprintf(HDR.FILE.stderr,'Error SCPOPEN: length=%i of channel %i different to length=%i of channel 1 \n',size(x,1),k,size(S2,1));
                                                return;
                                        else
	                                        fprintf(HDR.FILE.stderr,'Error SCPOPEN: Huffman decoding failed (%i) \n',size(x,1));
	    					HDR.data=S2;
						return;
                                        end;
                                end;
                                        
                        elseif HDR.SCP2.NHT~=19999,
                        	%% OBSOLETE %%
                                fprintf(HDR.FILE.stderr,'Error SOPEN SCP-ECG: user specified Huffman Table not supported\n');
                                HDR.SCP = SCP;
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
                        
                        if section.ID==5,
                                HDR.SCP5 = SCP;
                                HDR.SCP5.data = S2;
                                HDR.SampleRate = SCP.SampleRate;
                                
                        elseif section.ID==6,
                                HDR.SCP6 = SCP;
                                HDR.SampleRate = SCP.SampleRate;
                                HDR.PhysDim = repmat({HDR.SCP6.PhysDim},HDR.NS,1);
                                HDR.data = S2;

                                if HDR.FLAG.bimodal_compression,
                                	%% FIXME: THIS IS A HACK - DO NOT KNOW WHETHER IT IS CORRECT. 
                                %	HDR.FLAG.bimodal_compression = isfield(HDR,'SCP5') & isfield(HDR,'SCP4');
                                	HDR.FLAG.bimodal_compression = isfield(HDR,'SCP4');
				end; 
                                if HDR.FLAG.bimodal_compression,
					if isfield(HDR,'SCP5')
	                                        F = HDR.SCP5.SampleRate/HDR.SCP6.SampleRate;
        	                                HDR.SampleRate = HDR.SCP5.SampleRate;
                	                        HDR.FLAG.F = F;
					else
                	                        HDR.FLAG.F = 1;
					end;
                                        
                                        tmp=[HDR.SCP4.PA(:,1);HDR.LeadPos(1,2)]-[1;HDR.SCP4.PA(:,2)+1];
                                        if ~all(tmp==floor(tmp))
                                                tmp,
                                        end;
                                        t  = (1:HDR.N) / HDR.SampleRate;
                                        S1 = zeros(HDR.N, HDR.NS);
                                        
                                        p  = 1;
                                        k2 = 1;
                                        pa = [HDR.SCP4.PA;NaN,NaN];
                                        flag = 1;
                                        %% FIXME: accu undefined ##
					accu = 0;
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
                                                        S1(k1,:) = ((F-1)*accu + S2(fix(k2),:)) / F;
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
                                
                                if HDR.FLAG.ReferenceBeat & ~isfield(HDR,'SCP5') 
	                                fprintf(HDR.FILE.stderr,'Warning SOPEN SCP-ECG: Flag ReferenceBeat set, but no section 5 (containing the reference beat) is available\n');
                                elseif HDR.FLAG.ReferenceBeat,

                                	tmp_data = HDR.SCP5.data*(HDR.SCP5.Cal/HDR.SCP6.Cal); 
                                        for k = find(~HDR.SCP4.type(:,1)'),
                                                t1 = (HDR.SCP4.type(k,2):HDR.SCP4.type(k,4));
                                                t0 = t1 - HDR.SCP4.type(k,3) + HDR.SCP4.fc0;
                                                S2(t1,:) = S2(t1,:) + tmp_data(t0,:); 
                                        end;
                                end;
	                        HDR.data  = S2;
                        end;

                elseif section.ID==7, 
                        HDR.SCP7.byte1   = fread(fid,1,'uint8');    
                        HDR.SCP7.Nspikes = fread(fid,1,'uint8');    
                        HDR.SCP7.meanPPI = fread(fid,1,'uint16');    
                        HDR.SCP7.avePPI  = fread(fid,1,'uint16');    
                        
                        for k=1:HDR.SCP7.byte1,
                                HDR.SCP7.RefBeat{k} = fread(fid,16,'uint8');    
                                %HDR.SCP7.RefBeat1 = fread(fid,16,'uint8');    
                        end;
                        
                        for k=1:HDR.SCP7.Nspikes,
                                tmp = fread(fid,16,'uint16');    
                                tmp(1,2) = fread(fid,16,'int16');    
                                tmp(1,3) = fread(fid,16,'uint16');    
                                tmp(1,4) = fread(fid,16,'int16');    
                                HDR.SCP7.ST(k,:) = tmp;
                        end;
                        for k=1:HDR.SCP7.Nspikes,
                                tmp = fread(fid,6,'uint8');    
                                HDR.SCP7.ST2(k,:) = tmp;
                        end;
                        HDR.SCP7.Nqrs = fread(fid,1,'uint16');    
                        HDR.SCP7.beattype = fread(fid,HDR.SCP7.Nqrs,'uint8');    
                        
                        HDR.SCP7.VentricularRate = fread(fid,1,'uint16');    
                        HDR.SCP7.AterialRate = fread(fid,1,'uint16');    
                        HDR.SCP7.QTcorrected = fread(fid,1,'uint16');    
                        HDR.SCP7.TypeHRcorr = fread(fid,1,'uint8');    
                        
                        len = fread(fid,1,'uint16');
                        tag = 255*(len==0); 
                        k1 = 0;
                        while tag~=255,
                                tag = fread(fid,1,'uchar');    
                                len = fread(fid,1,'uint16');    
                                field = fread(fid,[1,len],'uchar');    
                                
                                if tag == 0,	
                                        HDR.Patient.LastName = char(field);
                                elseif tag == 1,
                                        
                                end;
                        end;
                        HDR.SCP7.P_onset = fread(fid,1,'uint16');    
                        HDR.SCP7.P_offset = fread(fid,1,'uint16');    
                        HDR.SCP7.QRS_onset = fread(fid,1,'uint16');    
                        HDR.SCP7.QRS_offset = fread(fid,1,'uint16');    
                        HDR.SCP7.T_offset = fread(fid,1,'uint16');    
                        HDR.SCP7.P_axis = fread(fid,1,'uint16');    
                        HDR.SCP7.QRS_axis = fread(fid,1,'uint16');    
                        HDR.SCP7.T_axis = fread(fid,1,'uint16');    
                        
                elseif section.ID==8, 
                        tmp = fread(fid,9,'uint8');    
                        HDR.SCP8.Report = tmp(1);    
                        HDR.SCP8.Time = [[1,256]*tmp(2:3),tmp(4:8)'];    
                        HDR.SCP8.N = tmp(9);    
                        for k = 1:HDR.SCP8.N,
                                ix  = fread(fid,1,'uint8');
                                len = fread(fid,1,'uint16');
                                tmp = fread(fid,[1,len],'uchar');    
                                HDR.SCP8.Statement{k,1} = char(tmp);    
                        end
                        
                %elseif section.ID==9, 
                %        HDR.SCP9.byte1 = fread(fid,1,'uint8');    
                        
                elseif section.ID==10, 
                        tmp = fread(fid,2,'uint16');
                        HDR.SCP10.NumberOfLeads = tmp(1);
                        HDR.SCP10.ManufacturerCode = tmp(2);
                        for k = []; 1:HDR.SCP10.NumberOfLeads,
                                tmp = fread(fid,2,'uint16')
                                LeadId = tmp(1); 
                                LeadLen = tmp(2); 
                                tmp = fread(fid,LeadLen/2,'uint16');    
                                HDR.SCP10.LeadId(k)=LeadId; 
                                HDR.SCP10.LeadLen(k)=LeadLen; 
                                HDR.SCP10.Measurements{k}=tmp; 
                        end;
                        
                elseif section.ID==11, 
                        bytes = fread(fid,11,'uint8');    
                        HDR.SCP11.T0 = [bytes(2)*256+bytes(3), bytes(4:8)];
                        HDR.SCP11.Confirmed = bytes(1);
                        HDR.SCP11.NumberOfStatements = bytes(9);
                        for k = 1:HDR.SCP11.NumberOfStatements,
                                SeqNo = fread(fid,1,'uint8');
                                len11 = fread(fid,1,'uint16');
                                typeID = fread(fid,1,'uint8');
                                Statement = fread(fid,len11-1,'uint8');
                                HDR.SCP11.Statement.SeqNo(k) = SeqNo;         
                                HDR.SCP11.Statement.len11(k) = len11;         
                                HDR.SCP11.Statement.typeID(k) = typeID;         
                                HDR.SCP11.Statement.Statement{k} = Statement;         
                        end;
                end;
                
		if ~section.Length,
			HDR.ERROR.status  = -1; 
			HDR.ERROR.message = 'Error SCPOPEN: \n';
			return;
		end;			
        end;

        HDR.SPR  = size(HDR.data,1);
        HDR.NRec = 1;
        HDR.AS.endpos = HDR.SPR;
        
        HDR.FILE.OPEN = 0; 
        HDR.FILE.POS  = 0;
        HDR.TYPE = 'native'; 
        fclose(HDR.FILE.FID);

else    % writing SCP file 

	NSections = 12;
        SectIdHdr = zeros(1,16); 
        VERSION = round(HDR.VERSION*10); 
        if ~any(VERSION==[10,13,20])
                fprintf(HDR.FILE.stderr,'Warning SCPOPEN(WRITE): unknown Version number %4.2f\n',HDR.VERSION);
              	VERSION = 20; 
        end;
	SectIdHdr(9:10) = VERSION; % Section and Protocol version number

        POS = 6; B = zeros(1,POS);
        for K = 0:NSections-1;
                b = [];
                if K==0,
                        % SECTION 0
                        b = [SectIdHdr(1:10),'SCPECG', zeros(1,NSections*10)];
	                b(16+(7:10)) = s4b(POS+1);

                elseif K==1,
                        % SECTION 1
                        b = SectIdHdr;
                        % tag(1),len(1:2),field(1:len)
                        if isfield(HDR.Patient,'Name'),
				b = [b, 0, s2b(length(HDR.Patient.Name)), HDR.Patient.Name];
			end;	
                        if isfield(HDR.Patient,'Id'),
				b = [b, 2, s2b(length(HDR.Patient.Id)), HDR.Patient.Id];
			end;
                        if isfield(HDR.Patient,'Age'),
				%b = [b, 4, s2b(3), s2b(HDR.Patient.Age),1];  %% use birthday instead
			end;
                        if isfield(HDR.Patient,'Birthday'),
	                        b = [b, 5, s2b(4), s2b(HDR.Patient.Birthday(1)),HDR.Patient.Birthday(2:3)];
			end;      
                        if isfield(HDR.Patient,'Height'),
                        if ~isnan(HDR.Patient.Height),
                        	b = [b, 6, s2b(3), s2b(HDR.Patient.Height),1];
                        end;
                        end;	
                        if isfield(HDR.Patient,'Weight'),
                        if ~isnan(HDR.Patient.Weight),
                        	b = [b, 7, s2b(3), s2b(HDR.Patient.Weight),1];
                        end;
                        end;	
	                if isfield(HDR.Patient,'Sex'),
	                if ~isempty(HDR.Patient.Sex)
                        	sex = HDR.Patient.Sex;
                        	if 0,
                        	elseif isnumeric(sex),           sex = sex(1); 
                        	elseif strncmpi(sex,'male',1);   sex = 1; 
                        	elseif strncmpi(sex,'female',1); sex = 2;
                        	else sex = 9; % unspecified
                        	end;
	                      	b = [b, 8, s2b(1), sex];
	                end;      	
                        end;	
                        if isfield(HDR.Patient,'Race'),
 	                       	b = [b, 9, s2b(1), HDR.Patient.Race(1)];
                        end;	
                        if isfield(HDR.Patient,'BloodPressure'),
                        	b = [b, 11, s2b(2), s2b(HDR.Patient.BloodPressure.Systolic)];
                        	b = [b, 12, s2b(2), s2b(HDR.Patient.BloodPressure.Diastolic)];
                        end;
                        
                        %% Tag 14
                        tag14.AnalyzingProgramRevisionNumber = ['',char(0)];
                        tag14.SerialNumberAcqDevice = ['',char(0)];
                        tag14.AcqDeviceSystemSoftware = ['',char(0)];
                        tag14.SCPImplementationSoftware = ['BioSig4OctMat v 1.76+',char(0)];	
                        tag14.ManufactureAcqDevice = ['',char(0)];	
                        t14 = [zeros(1,35), length(tag14.AnalyzingProgramRevisionNumber),tag14.AnalyzingProgramRevisionNumber,tag14.SerialNumberAcqDevice,tag14.AcqDeviceSystemSoftware,tag14.SCPImplementationSoftware,tag14.ManufactureAcqDevice];
                        t14(8)  = 255;  % Manufacturer
                        %%% ### FIXME ###  t14(9:14) = % cardiograph model
                        t14(15) = VERSION;        % Version
                        t14(16) = hex2dec('A0');  % Demographics and ECG rhythm data" (if we had also the reference beats we should change it in 0xC0).
                        t14(18) = hex2dec('D0');  % Capabilities of the ECG Device: 0xD0 (acquire, print and store). 
                        b   = [b, 14, s2b(length(t14)), t14];
                         
                        b = [b, 25, s2b(4), s2b(HDR.T0(1)),HDR.T0(2:3)];
                        b = [b, 26, s2b(3), HDR.T0(4:6)];
                        if ~any(isnan(HDR.Filter.HighPass))
	                        b = [b, 27, s2b(2), s2b(round(HDR.Filter.HighPass(1)*100))];
                        end; 
                        if ~any(isnan(HDR.Filter.LowPass))
	                	b = [b, 28, s2b(2), s2b(round(HDR.Filter.LowPass(1)))];
			end;
                        b = [b, 255, 0, 0];	% terminator
			b = b + (b<0)*256;

                elseif K==3,
                        % SECTION 3
                        b = [SectIdHdr,HDR.NS,4+HDR.NS*8];
                        if ~isfield(HDR,'LeadIdCode'), HDR.LeadIdCode = zeros(1,HDR.NS); end; 
                        if (numel(HDR.LeadIdCode)~=HDR.NS); warning('HDR.LeadIdCode does not have HDR.NS elements'); end; 
                        if any(HDR.LeadIdCode>255), warning('invalid LeadIdCode'); end;
                        for k = 1:HDR.NS,
                                b = [b, s4b(1), s4b(HDR.SPR*HDR.NRec), mod(HDR.LeadIdCode(k),256)];
                        end;

                elseif K==6,
                        % SECTION 6
                        Cal = full(HDR.Calib(2:end,:)); Cal = Cal - diag(diag(Cal)); 
                        if any(Cal(:))
                                fprintf(HDR.FILE.stderr,'Calibration is not a diagonal matrix.\n\tThis can result in incorrect scalings.\n'); 
                        end;
                        Cal = full(diag(HDR.Calib(2:end,:)));
                        if any(Cal~=Cal(1)), 
                                fprintf(HDR.FILE.stderr,'scaling information is not equal for all channels; \n\tThis is not supported by SCP and can result in incorrect scalings.\n'); 
                        end;
                        [tmp,scale1] = physicalunits(HDR.PhysDim{1});
                        [tmp,scale2] = physicalunits('nV');
                        b = [SectIdHdr, s2b(round(Cal(1)*scale1/scale2)), s2b(round(1e6/HDR.SampleRate)), 0, 0];
                        for k = 1:HDR.NS,
                                b = [b, s2b(HDR.SPR*HDR.NRec*2)];
                        end;
                        data = HDR.data(:);
                        data = data + (data<0)*2^16;
                        tmp  = s2b(round(data))';
                        b = [b,tmp(:)'];
                else
                        b = [];
                end;
                if (length(b)>0)
                        if mod(length(b),2), % align to multiple of 2-byte blocks
                                b = [b,0];
                        end;
                        if (length(b)<16), fprintf(HDR.FILE.stderr,'section header %i less then 16 bytes %i', K,length(b)); end;
                        b(3:4) = s2b(K);
                        b(5:8) = s4b(length(b));
                        %b(1:2)= s4b(crc);
                        b(1:2) = s2b(crc16eval(b(3:end)));
                        % section 0: startpos in pointer field 
	                B(22+K*10+(7:10)) = s4b(POS+1);
                end;
                B = [B(1:POS),b]; 
                % section 0 pointer field 
                B(22+K*10+(1:2)) = s2b(K);
                B(22+K*10+(3:6)) = s4b(length(b)); % length
                POS = POS + length(b);
        end
        B(3:6) = s4b(POS);      % length of file
        B(7:8) = s2b(crc16eval(B(9:22+NSections*10)));  % CRC of Section 0

        B(1:2) = s2b(crc16eval(B(3:end)));

        %        fwrite(fid,crc,'int16');
        count = fwrite(fid,B,'uchar');
        fclose(fid); 
end
end  %% scpopen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Auxillary functions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b2 = s2b(i)
	% converts 16bit into 2 bytes
	b2 = [bitand(i,255),bitand(bitshift(i,-8),255)];
	return; 
end	%%%%% s2b %%%%%


function b4 = s4b(i)
	% converts 32 bit into 4 bytes
	b4 = [s2b(bitand(i,2^16-1)),s2b(bitand(bitshift(i,-16),2^16-1)) ];
	return; 
end	%%%%% s4b %%%%%


function T = makeSubTree(T,bc,len,val)
	if (len==0)
		T.idxTable = val;
		return 
	end; 
if 1,
	b = bitand(bc,1)+1;
	if ~isfield(T,'branch') T.branch = {[],[]}; end;
	T.branch{b} = makeSubTree(T.branch{b},bitshift(bc,-1),len-1,val);
else 
	if bitand(bc,1)
		if ~isfield(T,'node1'),T.node1 = []; end;  
		T.node1 = makeSubTree(T.node1,bitshift(bc,-1),len-1,val);
	else
		if ~isfield(T,'node0'),T.node0 = []; end;  
		T.node0 = makeSubTree(T.node0,bitshift(bc,-1),len-1,val); 
	end; 
end,
	return; 
end	%%%%% makeSubTree %%%%%

function T = makeTree(HT)
	T = []; 
	for k1 = 1:size(HT,1)
		for k2 = 1:HT(k1,1) % CodeLength
			T = makeSubTree(T,HT(k1,5),HT(k1,1),k1); 
		end;  
	end; 
%save matlab T,pause
	return; 
end	%%%%% makeTree %%%%%

function outdata = DecodeHuffman(HTrees,HTs,indata,outlen)
	ActualTable = 1; 
	k1 = 1; r=0; 
	k2 = 0; 

	if ((outlen>0) && isfinite(outlen))
		outdata = repmat(NaN,outlen,1);
	else 	
		outdata = [0;0];  %% make it a column vector
	end; 	 
	Node    = HTrees{ActualTable}; 
	while ((k1*8+r <= 8*length(indata)) && (k2<outlen))
		if ~isfield(Node,'idxTable')
			r = r+1; if (r>8), k1=k1+1; r=1; end;
if 1,
			b = bitand(bitshift(indata(k1),r-8),1)+1;
			if ~isempty(Node.branch{b})
				Node = Node.branch{b};
			else 
				fprintf(2,'Warning SCPOPEN: empty node in Huffman table\n');
			end; 		 
else
			if bitand(bitshift(indata(k1),r-8),1)
				if isfield(Node,'node1') 
					Node = Node.node1;
				else 
					fprintf(2,'Warning SCPOPEN: empty node in Huffman table\n');
				end; 		 
			else 	 
				if isfield(Node,'node0') 
					Node = Node.node0; 
				else 
					fprintf(2,'Warning SCPOPEN: empty node in Huffman table\n');
				end; 
			end; 		
end;
		end; 

		if isfield(Node,'idxTable')
			TableEntry = HTs{ActualTable}(Node.idxTable,:);
			dlen = TableEntry(2)-TableEntry(1); 
			if (~TableEntry(3))
				ActualTable = TableEntry(4); 
			elseif (dlen~=0) 
				acc = 0;
				for k3 = 1:dlen,
					r = r+1; if (r>8), k1=k1+1; r=1; end;
					acc = 2*acc + bitand(bitshift(indata(k1),r-8), 1);
				end;
				if (acc>=bitshift(1,dlen-1))
					acc = acc - bitshift(1,dlen);  
				end; 	
				k2 = k2+1;
				outdata(k2) = acc; 
			else
				k2 = k2+1;
				outdata(k2)=TableEntry(4);
			end;
			Node = HTrees{ActualTable}; 
		end;
	end;
	return; 
end %%%%%%%% DecodeHuffman %%%%%%%%%

function crc16 = crc16eval(D)
% CRC16EVAL cyclic redundancy check with the polynomiaL x^16+x^12+x^5+1  
% i.e. CRC-CCITT http://en.wikipedia.org/wiki/Crc16 

	D = uint16(D);

	crchi = 255;
	crclo = 255;

	t = '00102030405060708191a1b1c1d1e1f112023222524272629383b3a3d3c3f3e32434041464744454a5b58595e5f5c5d53626160676665646b7a79787f7e7d7c74858687808182838c9d9e9f98999a9b95a4a7a6a1a0a3a2adbcbfbeb9b8bbbab6c7c4c5c2c3c0c1cedfdcdddadbd8d9d7e6e5e4e3e2e1e0effefdfcfbfaf9f8f9181b1a1d1c1f1e110003020504070608393a3b3c3d3e3f30212223242526272b5a59585f5e5d5c53424140474645444a7b78797e7f7c7d72636061666764656d9c9f9e99989b9a95848786818083828cbdbebfb8b9babbb4a5a6a7a0a1a2a3afdedddcdbdad9d8d7c6c5c4c3c2c1c0cefffcfdfafbf8f9f6e7e4e5e2e3e0e1e';
	crc16htab = hex2dec(reshape(t,2,length(t)/2)');

	t = '0021426384a5c6e708294a6b8cadceef31107352b594f7d639187b5abd9cffde62432001e6c7a4856a4b2809eecfac8d53721130d7f695b45b7a1938dffe9dbcc4e586a740610223cced8eaf48690a2bf5d4b79671503312fddcbf9e79583b1aa687e4c522036041ae8feccd2a0b684997b6d5f4133251709fbeddfc1b3a597888a9caeb0c2d4e6f80a1c2e304254667b998fbda3d1c7f5eb190f3d235147756eacba8896e4f2c0de2c3a08166472405dbfa99b85f7e1d3cd3f291b0577615344c6d0e2fc8e98aab44650627c0e182a37d5c3f1ef9d8bb9a75543716f1d0b3922e0f6c4daa8be8c926076445a283e0c11f3e5d7c9bbad9f81736557493b2d1f0';
	crc16ltab = hex2dec(reshape(t,2,length(t)/2)');

	for k = 1:length(D),
		ix = double(bitxor(crchi,D(k)))+1;
		crchi = bitxor(crclo,crc16htab(ix));
		crclo = crc16ltab(ix);
	end;
	crc16 = crchi*256+crclo;
	
end	%%%%% crc16eval %%%%%


