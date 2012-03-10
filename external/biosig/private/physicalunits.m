function [out,scale] = physicalunits(arg1)
% PHYSICALUNITS converts PhysDim inte PhysDimCode and vice versa
% according to Annex A of FEF Vital Signs Format [1]
%
%   HDR = physicalunits(HDR); 
%	adds HDR.PhysDim or HDR.PhysDimCode, if needed
%
%   PhysDim = physicalunits(PhysDimCode); 
%	converts Code of PhysicalUnits into descriptive physical units 
%
%   PhysDimCode = physicalunits(PhysDim); 
%	converts descriptive units into Code for physical units. 
%
%   [..., scale] = physicalunits(...);
%	scale contains the scaling factor of the decimal prefix 
%
% see also: SLOAD, SOPEN, doc/DecimalFactors.txt, doc/units.csv
%
% Reference(s): 
% [1] CEN/TC251/PT40 (2001)	
% 	File Exchange Format for Vital Signs - Annex A 
% [2] The Unified Code for Units of Measure (UCUM)
% 	http://aurora.regenstrief.org/UCUM/ucum.html
	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.

%	$Id$
%	Copyright (C) 2005,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/



global BIOSIG_GLOBAL;
if ~isfield(BIOSIG_GLOBAL,'ISLOADED')
	BIOSIG_GLOBAL.ISLOADED = 0 ; 
end; 
if ~BIOSIG_GLOBAL.ISLOADED; 
	f=which('getfiletype.m'); 	% identify path to biosig
	[p,f,e]=fileparts(f);
	[p,f,e]=fileparts(p);
	
	%%%---------- Load Decimal factors ------------%%%
	BIOSIG_GLOBAL.ISLOADED = 0 ; 
	[fid,msg]  = fopen(fullfile(p,'doc','DecimalFactors.txt'),'rt','ieee-le');
        line = fgetl(fid);
        N1   = 0; N2 = 0; 
        while ~feof(fid), 
		if ~strncmp(line,'#',1),
			N1 = N1 + 1;
               		[n,v,s] = str2double(line);
               		n = n(~v);
               		DecimalFactor.Code(N1,1) = n(2);
               		DecimalFactor.Cal(N1,1) = n(1);
               		s = s(~~v);
               		if any(v)
                      		DecimalFactor.Name(N1,1) = s(1);
                      		DecimalFactor.Prefix(N1,1) = s(2);
                       	end;
                end;	
                line = fgetl(fid);
        end;
	fclose(fid);
	BIOSIG_GLOBAL.DecimalFactor = DecimalFactor;

	%%%---------- Physical units ------------%%%
	fid = fopen(fullfile(p,'doc','units.csv'));
        line = fgetl(fid);
        N1 = 0; N2 = 0;
        while ~feof(fid),
		N2 = N2 + 1;
		if ~strncmp(line,'#',1),
			ix = mod(cumsum(line==char(34)),2); %% " 
			tmp = line; 
			tmp(~~ix) = ' '; 
			ix  = find(tmp==',');
			if (length(ix)~=3)
				fprintf(2,'Warning: line (%3i: %s) not valid\n',N2,line);
			else
				t1 = line(1:ix(1)-1);
				t2 = line(ix(1)+1:ix(2)-1);
				t3 = line(ix(2)+1:ix(3)-1);
				t4 = line(ix(3)+1:end);
				Code = str2double(t1);
				if ~isempty(Code)
					N1 = N1 + 1;
					UnitsOfMeasurement.Code(N1,1)   = Code;
					ix = min(find([t2, '[' ] == '['))-1;
					UnitsOfMeasurement.Symbol{N1,1} = char(t2(2:ix-1));
				end;	
              		end;
                end;	
                line = fgetl(fid);
        end;
	fclose(fid);
	BIOSIG_GLOBAL.Units = UnitsOfMeasurement;

	%%%---------- Load x73-UCUM Table ------------%%%
	fid=fopen('IEEEandUCUM.1b.txt','rt');
	if fid>1,
		r = char(fread(fid,[1,inf],'uint8')); 
		fclose(fid);
		while strncmp(r,'#',1)
			[line,r] = strtok(r,[10,13]); 
		end;
		K1 = 1; 
		[line,r] = strtok(r,[10,13]);
		while ~isempty(r),
			ix = [0,find(line == 9),length(line)+1]; % tab
			for K2 = 1:length(ix)-2, 
				if ix(K2)+2>ix(K2+1)-2,
					tab{K1,K2} = ' ';	
				else	
					tab{K1,K2} = line(ix(K2)+2:ix(K2+1)-2); 
				end; 	
			end;
			K2 = length(ix)-1; 
			tab{K1,6} = str2double(line(ix(K2)+1:ix(K2+1)-1));
			 
			[line,r] = strtok(r,[10,13]);
			K1 = K1+1; 
		end;
		BIOSIG_GLOBAL.x73_UCUM = tab;
	end; 	

	BIOSIG_GLOBAL.ISLOADED = 1;
end; 


if isstruct(arg1) 
	HDR = arg1;
	if 0, 
	elseif  isfield(HDR,'PhysDim') && isfield(HDR,'PhysDimCode')
		[PhysDimCode,scale] = physicalunits(HDR.PhysDim);
		if any(PhysDimCode(:)~=HDR.PhysDimCode(:))
			warning('PhysDim and PhysDimCode differ');
			PhysDim',HDR.PhysDim'
		end;
	elseif ~isfield(HDR,'PhysDim') && isfield(HDR,'PhysDimCode')
		[HDR.PhysDim, scale] = physicalunits(HDR.PhysDimCode);
	elseif  isfield(HDR,'PhysDim') % ~isfield(HDR,'PhysDimCode')
		[HDR.PhysDimCode,scale] = physicalunits(HDR.PhysDim);
	elseif ~isfield(HDR,'PhysDim') && ~isfield(HDR,'PhysDimCode')
		HDR.PhysDimCode = zeros(HDR.NS,1); 
		HDR.PhysDim = repmat({'?'},HDR.NS,1); 
		scale = ones(HDR.NS,1); 
		fprintf(HDR.FILE.stderr,'Neither PhysDim nor PhysDimCode defined in file %s\n',HDR.FileName);
	end;
	out = HDR;	

elseif isnumeric(arg1)
	s = mod(arg1,32);
	n = bitand(arg1,2^16-32);
	scale = repmat(NaN,size(arg1));
	PhysDim = repmat({'?'},size(arg1));
	for k = 1:length(n); 
	if any(BIOSIG_GLOBAL.DecimalFactor.Code==s(k)),
		t1 = BIOSIG_GLOBAL.DecimalFactor.Prefix{BIOSIG_GLOBAL.DecimalFactor.Code==s(k)};
		ix = find(BIOSIG_GLOBAL.Units.Code==n(k));
		if ~isempty(ix)
			t2 = BIOSIG_GLOBAL.Units.Symbol{ix};
			PhysDim{k} = [t1,t2];
			scale(k) = BIOSIG_GLOBAL.DecimalFactor.Cal(BIOSIG_GLOBAL.DecimalFactor.Code==s(k));
		end; 	
	end;
	end; 
	out = PhysDim;
	
elseif ischar(arg1) || iscell(arg1) 
	if iscell(arg1)
                N    = length(arg1); 
        elseif ischar(arg1)     
                N = size(arg1,1);
                arg1 = cellstr(arg1);
        end;
        j = 1:length(arg1);
        VER = version;

        % This does not work for Octave 2.9.9 and earlier versions. 
        % Uncomment the next line to fix it. It is slower but correct. 
	[arg1,i,j] = unique(arg1);

	Code = zeros(length(arg1),1); 	% default value is 0 (unknown)
	for k=1:length(arg1); 
		unit = deblank(arg1{k});
		if length(unit)>0,
			if (unit(1)==181), unit(1)='u'; end; 
		else 

		end;
                
                % this lookup table contains first the most widely used
                % physical units 
                if ~length(unit),
			Code(k) = 0;
		elseif strcmpi(unit,'-')	% dimensionless
               		Code(k) = 512;
		elseif strcmpi(unit,'1')	% dimensionless
               		Code(k) = 512;
		elseif strcmpi(unit,'percent')
               		Code(k) = 544;
		elseif strcmpi(unit,'%')
               		Code(k) = 544;
		elseif strcmpi(unit,'degree')
               		Code(k) = 736;
		elseif strcmpi(unit,'rad')
               		Code(k) = 768;
		elseif strcmp(unit,'Hz')
              		Code(k) = 2496;
		elseif strcmp(unit,'l/(min*m²)')
              		Code(k) = 2848;
		elseif strcmp(unit,'l/min')
              		Code(k) = 3072;
		elseif strcmp(unit,'mmHg')
              		Code(k) = 3872;
		elseif strcmp(unit,'dyne*s/cm^5')
              		Code(k) = 4128;
		elseif strcmp(unit,'V')
               		Code(k) = 4256;
		elseif strcmp(unit,'mV')
               		Code(k) = 4256 + 18;
		elseif strcmp(unit,'uV')
              		Code(k) = 4256 + 19;
		elseif strcmp(unit,'Ohm')
               		Code(k) = 4288;
		elseif strcmp(unit,'K')
              		Code(k) = 4384;
		elseif strcmp(unit,'°F')
              		Code(k) = 4416;
		elseif strcmp(unit,'°C')
              		Code(k) = 6048;
		elseif strcmp(unit,'m/s²')
              		Code(k) = 6624
		elseif strcmp(unit,'dyne*s*m²/cm^5')
              		Code(k) = 65440;
		elseif strcmp(unit,'l/m²')
              		Code(k) = 65472;
		elseif strcmp(unit,'ml/m²')
              		Code(k) = 65472+18;
                else 
                        % this is the general method to determine the
                        % Physical Dimension's code 
	       		ix = [];
			for k1=1:length(BIOSIG_GLOBAL.DecimalFactor.Code)
			for k2=1:length(BIOSIG_GLOBAL.Units.Code)
			if strcmp(unit,[BIOSIG_GLOBAL.DecimalFactor.Prefix{k1},BIOSIG_GLOBAL.Units.Symbol{k2}])
				ix = [ix,BIOSIG_GLOBAL.Units.Code(k2) + BIOSIG_GLOBAL.DecimalFactor.Code(k1)];
			end;
			end;
			end;

			%%%%% usis x73-UCUM table, col 3 and 4 %%%%
			ix2 = []; 
			ix3 = []; 
			if isfield(BIOSIG_GLOBAL,'x73_UCUM')
			for k1=1:length(BIOSIG_GLOBAL.DecimalFactor.Code)
			for k2=1:size(BIOSIG_GLOBAL.x73_UCUM,1)
			if strcmpi(unit,[BIOSIG_GLOBAL.DecimalFactor.Prefix{k1},BIOSIG_GLOBAL.x73_UCUM{k2,3}])
				ix2 = [ix2,BIOSIG_GLOBAL.x73_UCUM{k2,6} + BIOSIG_GLOBAL.DecimalFactor.Code(k1)];
			end;
			if strcmpi(unit,[BIOSIG_GLOBAL.DecimalFactor.Prefix{k1},BIOSIG_GLOBAL.x73_UCUM{k2,4}])
				ix3 = [ix3,BIOSIG_GLOBAL.x73_UCUM{k2,6} + BIOSIG_GLOBAL.DecimalFactor.Code(k1)];
			end;
			end;
			end;
			end; 

			if length(ix)==1,
				Code(k) = ix; 
			elseif length(ix)>1,
				ix = ix(~mod(ix,32));	% select those with no PreFix (code offset = 0) 
				if length(ix)==1,
					Code(k) = ix; 
				elseif sum(ix<64000)==1,
					Code(k) = ix(ix<64000); 
				else  
					warning('ambigous physical unit')
				end; 	
			elseif length(ix2)==1,
				Code(k) = ix2; 
			elseif length(ix3)==1,
				Code(k) = ix3; 
			end;	
                end;
                
        end;        
        Code  = Code(j); 
	scale = repmat(NaN,size(Code));
	for k = 1:numel(Code); 
		scale(k) = BIOSIG_GLOBAL.DecimalFactor.Cal(BIOSIG_GLOBAL.DecimalFactor.Code==bitand(Code(k),31));
	end;
        out = Code; 
        
end;
