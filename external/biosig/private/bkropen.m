function [BKR,s]=bkropen(arg1,arg3,arg4,arg5,arg6)
% BKROPEN opens BKR file
% However, it is recommended to use SOPEN instead .
% For loading whole data files, use SLOAD. 
%
% see also: SOPEN, SREAD, SSEEK, STELL, SCLOSE, SWRITE, SEOF

%	$Id$
%	Copyright (c) 1997-2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.

if isstruct(arg1),
	BKR=arg1;
	FILENAME=BKR.FileName;
	if BKR.FILE.OPEN,
		fseek(BKR.FILE.FID,0,'bof');	
	else
		BKR.FILE.FID = fopen(FILENAME,BKR.FILE.PERMISSION,'ieee-le');          
	end;
else
	FILENAME=arg1;
	BKR.FILE.FID = fopen(FILENAME,BKR.FILE.PERMISSION,'ieee-le');          
        BKR.FileName  = FILENAME;
        [pfad,file,FileExt] = fileparts(BKR.FileName);
        BKR.FILE.Name = file;
        BKR.FILE.Path = pfad;
        BKR.FILE.Ext  = FileExt(2:length(FileExt));
end;
if BKR.FILE.FID<0,
        fprintf(2,'Error BKROPEN: file %s not found.\n',FILENAME); 
        return;
end;
if ftell(BKR.FILE.FID)~=0,	% 
        fprintf(2,'Error: Fileposition is not 0\n');        	        
end;

BOOL ='int16';
ULONG='uint32'; 
FLOAT='float32';

if any(BKR.FILE.PERMISSION=='r'),
	BKR.FILE.OPEN = 1;
	fid = BKR.FILE.FID;
        %%%%% READ HEADER

        %VARIABLE	TYPE			#bytes	OFFSET	COMMENT
	BKR.VERSION     = fread(fid,1,'uint16');	%	2 Byte	0	Versionsnummer 
	if ((BKR.VERSION<=200) | (BKR.VERSION>207)) fprintf(2,'LOADBKR: WARNING  Version BKR Format %i',BKR.VERSION); end;
	BKR.NS          = fread(fid,1,'uint16');	%	2 Byte	2	Anzahl der Kanäle
	BKR.SampleRate  = fread(fid,1,'uint16');	%	2 Byte	4	Abtastfrequenz	
	BKR.NRec        = fread(fid,1,'uint32');	%	4 Byte	6	Anzahl der Trials
	BKR.SPR         = fread(fid,1,'uint32');	%	4 Byte	10	Anzahl Samples per Trial
	BKR.PhysMax(1:BKR.NS)     = fread(fid,1,'uint16');	%	2 Byte	14	Kalibrierspannung	
	BKR.DigMax(1:BKR.NS)      = fread(fid,1,'uint16');	%	2 Byte	16	Kalibrierwert
	Label           = fread(fid,[1,4],'uint8');	%	4 Byte	18	Elektrodencode
	BKR.Label 	= repmat({char(Label)},BKR.NS,1); 
	BKR.Filter.LowPass      = fread(fid,1, FLOAT);	%	4 Byte	22	untere Eckfrequenz
	BKR.Filter.HighPass     = fread(fid,1, FLOAT);	%	4 Byte	26	obere Eckfrequenz
	BKR.BKR.sref=fread(fid,1, ULONG);	%	4 Byte	30	Startzeitpunkt Referenz in Samples
	BKR.BKR.eref=fread(fid,1, ULONG);	%	4 Byte	34	Länge Referenz in Samples	
	BKR.BKR.sact=fread(fid,1, ULONG);	%	4 Byte	38	Startzeitpunkt Aktion in Samples
	BKR.BKR.eact=fread(fid,1, ULONG);	%	4 Byte	42	Länge Aktion in Samples
	BKR.FLAG.TRIGGERED      = fread(fid,1,BOOL);	%	2 Byte	46	flag für Trigger
	BKR.BKR.pre=fread(fid,1, ULONG);	%	4 Byte	48	Anzahl der Sampels vor dem Trigger
	BKR.BKR.pst=fread(fid,1, ULONG);	%	4 Byte	52	Anzahl der Sampels nach dem Trigger
	BKR.BKR.hav=fread(fid,1,BOOL);		%	2 Byte	56	flag für "horizontale" Mittelung
	BKR.BKR.nah=fread(fid,1, ULONG);	%	4 Byte	58	Anzahl der gemittelten Trials
	BKR.BKR.vav=fread(fid,1,BOOL);		%	2 Byte	62	flag für "vertikale" Mittelung	
	BKR.BKR.nav=fread(fid,1,'uint16');	%	2 Byte	64	Anzahl der gemittelten Kanäle
	BKR.BKR.cav=fread(fid,1,BOOL);		%	2 Byte	66	flag für Datenkomprimierung
	BKR.BKR.nac=fread(fid,1, ULONG);	%	4 Byte	68	Anzahl der gemittelten Samples
	BKR.FLAG.ref=fread(fid,4,BOOL);		%	2 Byte	72	flag: Common Average Reference
	%       loc=fread(fid,1,BOOL);		%	2 Byte	74	flag: Local Average Reference
	%       lap=fread(fid,1,BOOL);		%	2 Byte	76	flag: Laplace Berechnung
	%       wgt=fread(fid,1,BOOL);		%	2 Byte	78	flag: Weighted Average Reference
	BKR.BKR.pwr=fread(fid,1,BOOL);		%	2 Byte	80	flag: Leistung
	BKR.BKR.avr=fread(fid,1,BOOL);		%	2 Byte	82	flag: Mittelwert
	BKR.BKR.std=fread(fid,1,BOOL);		%	2 Byte	84	flag: Streuung
	BKR.BKR.bps=fread(fid,1,BOOL);		%	2 Byte	86	flag: Bandpaß
	BKR.BKR.erd=fread(fid,1,BOOL);		%	2 Byte	88	flag: ERD
	BKR.BKR.sig=fread(fid,1,BOOL);		%	2 Byte	90	flag: Signifikanz 
	BKR.BKR.coh=fread(fid,1,BOOL);		%	2 Byte	92	flag: Kohärenz
	BKR.BKR.spc=fread(fid,1,BOOL);		%	2 Byte	94	flag: Spectrum
	BKR.BKR.conf=fread(fid,1, FLOAT);	%	4 Byte	96	Konfidenz
	BKR.BKR.csp=fread(fid,1,BOOL);		%	2 Byte	100	flag: Kohärenz Leistungsspektrum
	BKR.BKR.erc=fread(fid,1,BOOL);		%	2 Byte	102	flag: ERC
	BKR.BKR.ham=fread(fid,1,BOOL);		%	2 Byte	104	flag: Hanning smoothed
	BKR.BKR.ann=fread(fid,1,BOOL);		%	2 Byte	106	flag: art. Neuronal. NW. Filter (ANN)
	niu=fread(fid,1,'uint16');	%	2 Byte 	108	ANN: Anzahl input units 
	nhu=fread(fid,1,'uint16');	%	2 Byte	110	ANN: Anzahl hidden units
	nlc=fread(fid,1, ULONG);	%	4 Byte	112	ANN: Anzahl Lernzyklen
	reg=fread(fid,1, FLOAT);	%	4 Byte	116	ANN: regression
	lco=fread(fid,1, FLOAT);	%	4 Byte 	120	ANN: Lernkoeffizient
	epo=fread(fid,1,'uint16');	%	2 Byte	124	ANN: Epoche
	BKR.BKR.rel=fread(fid,1,BOOL);		%	2 Byte	126	flag: ERC in Relativwerten
	wnd=fread(fid,1,'uint16');	%	2 Byte	128	Fenstertyp
	BKR.BKR.kal=fread(fid,1,BOOL);		%	2 Byte	130	flag: Kalman gefilterte Daten
	BKR.BKR.cwt=fread(fid,1,BOOL);		%	2 Byte	132	flag: kont. Wavelet transformtiert
	cwt_fmin=fread(fid,1, FLOAT);	%	4 Byte	134	unterste Frequenz, kont. Wavelettransform.
	cwt_fmax=fread(fid,1, FLOAT);	%	4 Byte	138	oberste Frequenz, kont. Wavelettransform.
	scales=fread(fid,1,'uint16');	%	2 Byte	142	Anzahl Frequenzbänder für kont. WT
	cwt_fe=fread(fid,1, FLOAT);	%	4 Byte	144	frequ. für Dt = Df = 1/2Öp
	cwt_start=fread(fid,1, ULONG);	%	4 Byte	148	Startsample für kont. WT Berechnung
        %-- NULL --	-------------	--------	152	-- Offset bis 512 Byte --
        if 1, 
                fread(fid,1024-152,'uint8');
        else
                fread(fid,512-152,'uint8');
                for i=1:BKR.NS,
                        eletyp(i)=fread(fid,1,'uchar');	%	1 Byte	512	Elektrode 1: Signalart (z.B: EEG)
                        elenum(i)=fread(fid,1,'uchar');	%	1 Byte	513	Elektrode 1: Kanalnr. für gleiche Signalart
                        ref(i)=fread(fid,1, FLOAT);	%	4 Byte	514	Referenzwert für Kanal 1
                end;
                if BKR.NS>85,
                        fprintf(2,'Warning BKRLOAD: Number of channels larger than 85; Header does not support more\n');
                end;
                fseek(fid,512-(BKR.NS*6),'cof');
        end;
        BKR.HeadLen = ftell(fid);
	if BKR.HeadLen~=1024,
	        fprintf(2,'Warning BKRLOAD: Length of Header does not confirm BKR-specification\n');
	end;
	%BKR.HeadLen = 1024;

	%%%%% Generate BKR-Struct according to biosig/doc/header.txt
	BKR.Dur=1/BKR.SampleRate;
	BKR.DigMin = -BKR.DigMax;
	BKR.PhysMin= -BKR.PhysMax;
	BKR.Cal=BKR.PhysMax./BKR.DigMax;
	BKR.Off=zeros(BKR.NS,1);
        BKR.Calib = sparse(2:BKR.NS+1,1:BKR.NS,BKR.Cal,BKR.NS+1,BKR.NS);
        %BKR.PhysDim = repmat({'µV'},BKR.NS,1);
        BKR.PhysDimCode = repmat(4275,BKR.NS,1);  % uV 
	tmp=sprintf('LowPass %4.1f Hz; HighPass %4.1f Hz; Notch ?',BKR.Filter.LowPass,BKR.Filter.HighPass);
	BKR.PreFilt=tmp; %ones(BKR.NS,1)*[tmp 32+zeros(1,80-length(tmp))];
	BKR.Filter.Notch    = nan; %h.notchfilter;

	BKR.AS.startrec = 0;
	BKR.AS.numrec = 0;
	BKR.AS.bpb = BKR.NS*2;	% Bytes per Block
	BKR.AS.spb = BKR.NS;	% Samples per Block
	BKR.FILE.POS = 0;

	if ~BKR.FLAG.TRIGGERED & (BKR.NRec>1);
		fprintf(BKR.FILE.stderr,'Warning: TriggerFlag in file %s was not set.\n',BKR.FileName);
		BKR.FLAG.TRIGGERED = 1;
	end;	
        BKR.FLAG.REFERENCE = 'unknown';
        if BKR.FLAG.ref(1), BKR.FLAG.REFERENCE = 'COM'; end;
        if BKR.FLAG.ref(2), BKR.FLAG.REFERENCE = 'LOC'; end;
        if BKR.FLAG.ref(3), BKR.FLAG.REFERENCE = 'LAP'; end;
        if BKR.FLAG.ref(4), BKR.FLAG.REFERENCE = 'WGT'; end;
                
        % THRESHOLD for Overflow detection
        BKR.AS.endpos = (BKR.FILE.size-BKR.HeadLen)/BKR.AS.bpb;

        BKR.data = fread(fid,[BKR.NS,inf],'int16')';
	fclose(fid);
        BKR.TYPE = 'native'; 

        % check whether Headerinfo fits to file length.
	if (BKR.FILE.size-BKR.HeadLen)~=BKR.SPR*BKR.NRec*BKR.NS*2,
		%[BKR.FILE.size,BKR.HeadLen,BKR.SPR,BKR.NRec,BKR.NS],
		%[BKR.FILE.size-BKR.HeadLen-BKR.SPR*BKR.NRec*BKR.NS*2],
		fprintf(2,'Warning BKROPEN: Header information in %s corrupted;',BKR.FileName);
		if BKR.NRec==1,
        		fprintf(2,'Data could be reconstructed.\n',BKR.FileName);
			BKR.SPR=(BKR.FILE.size-BKR.HeadLen)/(BKR.NRec*BKR.NS*2);
                elseif BKR.NRec==0,
        		fprintf(2,'Data could be reconstructed.\n',BKR.FileName);
			BKR.NRec=(BKR.FILE.size-BKR.HeadLen)/(BKR.SPR*BKR.NS*2);
                end;
        	if (BKR.FILE.size-BKR.HeadLen)~=BKR.SPR*BKR.NRec*BKR.NS*2,
        		fprintf(2,'Unable to reconstruct data.\n',BKR.FileName);
        	end;
        end;

        % look for Classlabel information
        if ~isfield(BKR,'Classlabel'),
                BKR.Classlabel = [];
        end;
        if ~isfield(BKR,'TRIG'),
                BKR.TRIG = [];
        end;
        tmp=fullfile(BKR.FILE.Path,[BKR.FILE.Name,'.mat']);
        if ~exist(tmp,'file'),
                tmp=fullfile(BKR.FILE.Path,[BKR.FILE.Name,'.MAT']);
        end
        x = [];
        if exist(tmp,'file'),
                x = load('-mat',tmp);
        end;
        if isfield(x,'header'),
                BKR.MAT  = x.header;
                if isfield(x.header,'Setup'), 
                        if isfield(x.header.Setup,'Bits'), 
                                BKR.Bits = x.header.Setup.Bits;
                                [datatyp, limits, datatypes] = gdfdatatype(BKR.Bits+255);
                                % THRESHOLD for Overflow detection
                                if ~isfield(BKR,'THRESHOLD')
	                                BKR.THRESHOLD = repmat(limits, BKR.NS, 1);
	                        end;         
                        end;
                end;
                if isfield(x.header,'Result') & isfield(x.header.Result,'Classlabel'),
                        BKR.Classlabel = x.header.Result.Classlabel;
                end;
                if isfield(x.header,'Paradigm')
                        if isempty(BKR.Classlabel) & isfield(x.header.Paradigm,'Classlabel')
                                BKR.Classlabel = x.header.Paradigm.Classlabel;
                        end;
                        BKR.BCI.Paradigm = x.header.Paradigm;
                        if isfield(BKR.BCI.Paradigm,'TriggerOnset');
                                BKR.TriggerOffset = BKR.BCI.Paradigm.TriggerOnset;
                        elseif isfield(BKR.BCI.Paradigm,'TriggerTiming');
                            %    BKR.BCI.Paradigm.TriggerTiming,
                                BKR.TriggerOffset = BKR.BCI.Paradigm.TriggerTiming;
                                fprintf(2,'Warning BKROPEN: Paradigm.TriggerOnset is unknown. Paradigm.TriggerTiming= %f ms is used instead\n',BKR.TriggerOffset);
                        end;
                end;

                if isfield(x.header,'PhysioRec'), % R. Leeb's data 
                        BKR.Label = x.header.PhysioRec;
                end;
                if isfield(x.header,'BKRHeader'), % R. Scherer Data 
                        if isfield(x.header.BKRHeader,'TO'),
                                BKR.T0 = x.header.BKRHeader.TO;
                        end;
                        if isfield(x.header.BKRHeader,'Label'),
                                BKR.Label = x.header.BKRHeader.Label;
                                ns = BKR.NS-size(BKR.Label,1);
                                if ns == 1;
                                        BKR.Label = strvcat(BKR.Label,'TRIGGER');
                                elseif ns > 1;
                                        BKR.Label = strvcat(BKR.Label,char(repmat('n.a.',ns,1)));
                                end;
                        end;
                end;
                if isfield(x.header,'Model'), % More 
                        if isfield(x.header.Model,'AnalogInput'), 
                                for k = 1:length(x.header.Model.AnalogInput),
                                        BKR.Filter.HighPass(k) = x.header.Model.AnalogInput{k}{5};
                                        BKR.Filter.LowPass(k)  = x.header.Model.AnalogInput{k}{6};
                                        BKR.Filter.Notch(k)    = strcmpi(x.header.Model.AnalogInput{k}{7},'on');

                                        BKR.MAT.Cal(k) = x.header.Model.AnalogInput{k}{3};
                                end
                        end;
                end;
                if ~isempty(strmatch('TRIGGER',BKR.Label))
                        BKR.AS.TRIGCHAN = BKR.NS; %strmatch('TRIGGER',H.Label); 
                end;
        end;
        if 1; %~isfield(BKR,'Classlabel'),
                tmp=fullfile(BKR.FILE.Path,[BKR.FILE.Name,'.par']);
                if ~exist(tmp,'file'),
                	tmp=fullfile(BKR.FILE.Path,[BKR.FILE.Name,'.PAR']);
                end
                if exist(tmp,'file'),
                        BKR.Classlabel = load(tmp);
                end;
        end;

        %%% Artifact Selection files 
        tmp1=fullfile(BKR.FILE.Path,[BKR.FILE.Name,'.sel']);
        if ~exist(tmp1,'file'),
                tmp1=fullfile(BKR.FILE.Path,[BKR.FILE.Name,'.SEL']);
        end
        tmp2 = fullfile(BKR.FILE.Path,[BKR.FILE.Name,'_artifact.mat']);
        SW   = (exist(tmp1,'file')>0) + 2*(exist(tmp2,'file')>0);
        if SW == 0, 
        elseif SW == 1,
                if exist('OCTAVE_VERSION')>5
                        BKR.ArtifactSelection = load('-ascii',tmp1);
                else
                        BKR.ArtifactSelection = load(tmp1);
                end;
        elseif SW == 2,
                if exist('OCTAVE_VERSION')>5
                        tmp = load('-mat',tmp2);
                else
                        tmp = load(tmp2);
                end;
                BKR.ArtifactSelection = tmp.artifact(:);
        elseif SW == 3,
                fprintf(BKR.FILE.stderr,'Warning BKROPEN: more than one ArtifactSelection files. File %s is used.\n',tmp1);
                if exist('OCTAVE_VERSION')>5
                        BKR.ArtifactSelection = load('-ascii',tmp1);
                else
                        BKR.ArtifactSelection = load(tmp1);
                end;
        end;
        if isfield(BKR,'ArtifactSelection'),
                if any(BKR.ArtifactSelection>1) | (length(BKR.ArtifactSelection)<length(BKR.Classlabel))
                        sel = zeros(size(BKR.Classlabel));
                        sel(BKR.ArtifactSelection) = 1;
                        BKR.ArtifactSelection = sel(:);
                end;
                BKR.ArtifactSelection = BKR.ArtifactSelection(:);
        end;
        
        if isfield(BKR.AS,'TRIGCHAN') % & isempty(BKR.EVENT.POS)
                if BKR.AS.TRIGCHAN <= BKR.NS, %size(BKR.data,2),
                        BKR.THRESHOLD(BKR.AS.TRIGCHAN,1:2) = [-1-2^15,2^15]; % do not apply overflow detection for Trigger channel 
                        TRIGon = gettrigger(double(BKR.data(:,BKR.AS.TRIGCHAN)));
                        %TRIGoff = gettrigger(-double(BKR.data(:,BKR.AS.TRIGCHAN)));
                        if isfield(BKR,'TriggerOffset')
                                TRIGon  = TRIGon - round(BKR.TriggerOffset/1000*BKR.SampleRate);
                        %        TRIGoff = TRIGoff - round(BKR.TriggerOffset/1000*BKR.SampleRate);
                        end;
                end;
                BKR.TRIG = TRIGon(:);
                BKR.EVENT.POS = TRIGon(:); %[TRIGon(:); TRIGoff(:)]; 
                BKR.EVENT.TYP = repmat(hex2dec('0300'),numel(TRIGon),1); %repmat(hex2dec('8300'),numel(TRIGoff),1)];
        end;
        if length(BKR.TRIG)~=length(BKR.Classlabel),
                % hack to deal with BCI22 data
                fprintf(2,'Warning BKROPEN: Number of triggers (%i) and number of Classlabels (%i) do not fit\n',length(BKR.TRIG),length(BKR.Classlabel));
                BKR.TRIG = [];
                BKR.Classlabel = [];
                BKR.ArtifactSelection = [];
        end;

        
elseif any(BKR.FILE.PERMISSION=='w'),
        
        BKR.FILE.OPEN = 2;		
        BKR.VERSION   = 207;
        BKR.TYPE      = 'BKR';
        if ~strcmpi(BKR.FILE.Ext,'BKR'),
                fprintf(2,'Warning BKROPEN-WRITE: file extionsion is not BKR.\n');
        end;
        if ~isfield(BKR,'SampleRate'),
                fprintf(2,'Error BKROPEN-WRITE: BKR.SampleRate not defined.\n');
                return;
        end;
        if ~isfield(BKR,'NS'),
                BKR.NS = 0; 	% unknown channel number ...
        end;
        if ~isfield(BKR,'SPR'),
                BKR.SPR = 0; 	% Unknown - Value will be fixed when file is closed. 
        end;
        if isfield(BKR,'NRec'),
                BKR.FLAG.TRIGGERED = BKR.NRec>1;
        else
                BKR.NRec = -1; 	% Unknown - Value will be fixed when file is closed. 
        end;
        if ~isfield(BKR,'PhysMax'), BKR.PhysMax = NaN; end;
        if isempty(BKR.PhysMax),    BKR.PhysMax = NaN; end;
        if ~isfield(BKR,'DigMax'),  BKR.DigMax  = NaN; end;
        if isnan(BKR.DigMax) | isempty(BKR.DigMax),
                BKR.DigMax = 2^15-1;   
        end;
        
        if any([BKR.NS==0,BKR.SPR==0,BKR.NRec<0,isnan([BKR.NRec,BKR.NS,BKR.SPR,BKR.DigMax,BKR.PhysMax,BKR.SampleRate])]), 	% if any unknown, ...	
                BKR.FILE.OPEN = 3;			%	... fix header when file is closed. 
        end;
        if ~isfield(BKR,'FLAG'),
                BKR.FLAG.UCAL = 0; 
        elseif ~isfield(BKR.FLAG,'UCAL'),
                BKR.FLAG.UCAL = 0; 
        end;

        tmp = round(BKR.PhysMax);
	fprintf(1,'Scaling error in file %s due to rounding of PhysMax is in the range of %f%%.\n',BKR.FileName, abs((BKR.PhysMax-tmp)/tmp)*100);
	BKR.PhysMax = tmp;

	count=fwrite(BKR.FILE.FID,BKR.VERSION,'short');	        % version number of header
	count=fwrite(BKR.FILE.FID,BKR.NS,'short');	        % number of channels
	count=fwrite(BKR.FILE.FID,BKR.SampleRate,'short');      % sampling rate
	count=fwrite(BKR.FILE.FID,BKR.NRec,'int32');            % number of trials: 1 for untriggered data
	count=fwrite(BKR.FILE.FID,BKR.SPR,'uint32');            % samples/trial/channel
	count=fwrite(BKR.FILE.FID,BKR.PhysMax,'short');		% Kalibrierspannung
	count=fwrite(BKR.FILE.FID,BKR.DigMax, 'short');		% Kalibrierwert
        
	count=fwrite(BKR.FILE.FID,zeros(4,1),'uint8');        
	if isfield(BKR,'Filter'),
		if ~isfield(BKR.Filter,'LowPass'),
			BKR.Filter.LowPass = NaN; 
		end;        
		if ~isfield(BKR.Filter,'HighPass'),
			BKR.Filter.HighPass= NaN; 
		end;
	else
		BKR.Filter.LowPass =NaN; 
		BKR.Filter.HighPass=NaN; 
        end;
	if all(BKR.Filter.LowPass(1)==BKR.Filter.LowPass)
		BKR.Filter.LowPass = BKR.Filter.LowPass(1);
	else
		BKR.Filter.LowPass = NaN;
	end;
	if all(BKR.Filter.HighPass(1)==BKR.Filter.HighPass)
		BKR.Filter.HighPass = BKR.Filter.HighPass(1);
	else
		BKR.Filter.HighPass = NaN;
	end;
	count=fwrite(BKR.FILE.FID,[BKR.Filter.LowPass,BKR.Filter.HighPass],'float'); 

	count=fwrite(BKR.FILE.FID,zeros(16,1),'uint8');         	% offset 30
	count=fwrite(BKR.FILE.FID,BKR.FLAG.TRIGGERED,'int16');	% offset 32
	count=fwrite(BKR.FILE.FID,zeros(24,1),'uint8');         	% offset 46
        
        if ~isfield(BKR,'FLAG')
                BKR.FLAG.REFERENCE='';
        end;
        if ~isfield(BKR.FLAG,'REFERENCE')
                BKR.FLAG.REFERENCE='';
        end;
        
        tmp  = [strcmp(BKR.FLAG.REFERENCE,'COM')|strcmp(BKR.FLAG.REFERENCE,'CAR'), strcmp(BKR.FLAG.REFERENCE,'LOC')|strcmp(BKR.FLAG.REFERENCE,'LAR'), strcmp(BKR.FLAG.REFERENCE,'LAP'), strcmp(BKR.FLAG.REFERENCE,'WGT')];
        
        fwrite(BKR.FILE.FID,tmp,BOOL); 		% offset 72 + 4*BOOL
                
	%speichert den rest des BKR-headers
	count = fwrite(BKR.FILE.FID,zeros(1024-80,1),'uint8');
	BKR.HeadLen = ftell(BKR.FILE.FID);
	if BKR.HeadLen~=1024,
		fprintf(2,'Error BKROPEN WRITE: HeaderLength is not 1024 but %i\n',BKR.HeadLen);
	end;
	BKR.FILE.POS  = 0;
	BKR.AS.endpos = 0;
	BKR.AS.bpb = BKR.NS*2;	% Bytes per Block
        BKR.AS.spb = BKR.NS;	% Samples per Block
        
        if isfield(BKR,'Classlabel');
                fid = fopen(fullfile(BKR.FILE.Path,[BKR.FILE.Name,'.par']),'w+b');  % force binary mode
                fprintf(fid,'%i\r\n',BKR.Classlabel);  % explicit 0x0d and 0x0a
                fclose(fid);
        end;
        if isfield(BKR,'ArtifactSelection');
                fid = fopen(fullfile(BKR.FILE.Path,[BKR.FILE.Name,'.sel']),'w+b');  % force binary mode
                fprintf(fid,'%i\r\n',BKR.ArtifactSelection);  % explicit 0x0d and 0x0a
                fclose(fid);
        end;
end;


