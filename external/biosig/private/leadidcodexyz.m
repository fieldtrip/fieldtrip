function [HDR] = leadidcodexyz(arg1)
% LeadIdCodeXYZ uses the Label information for computing the 
%       LeadIdCode and the XYZ position of the EEG Electrodes
%       according to Annex A of FEF Vital Signs Format [1]
%
%   HDR = leadidcodexyz(HDR); 
%	adds HDR.LeadIdCode and HDR.ELEC.XYZ, if needed. 
%
% see also: SLOAD, SOPEN, PHYSICALUNITS, doc/leadidtable_scpecg.txt, doc/elecpos.txt
%
% Reference(s): 
% [1] CEN/TC251/PT40 (2001)	
% 	File Exchange Format for Vital Signs - Annex A 
%
% Birbaumer, N. (2006). Brain-computer-interface research: Coming of age. Clinical Neurophysiology, 117:479–83. 
%  http://www.acns.org/pdfs/ACFDD46.pdf. 
% ACNS (2006). Guidelines for standard electrode position nomenclature. American Clinical
% Neurophysiology Society. http://www.acns.org/pdfs/ACFDD46.pdf.

% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.

%	$Id$
%	Copyright (C) 2006,2007,2008,2009 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


global BIOSIG_GLOBAL;
% BIOSIG_GLOBAL=[]; %%% used for debugging, only. 

if ~isfield(BIOSIG_GLOBAL,'Phi'); 
	BIOSIG_GLOBAL.ISLOADED_XYZ = 0 ; 
end; 
if ~isfield(BIOSIG_GLOBAL,'ISLOADED_XYZ')
	BIOSIG_GLOBAL.ISLOADED_XYZ = 0 ; 
end; 
if ~BIOSIG_GLOBAL.ISLOADED_XYZ; 
	f = which('getfiletype.m'); 	% identify path to biosig
	[p,f,e] = fileparts(f); 
	[p,f,e] = fileparts(p); 
	
        BIOSIG_GLOBAL.ISLOADED_XYZ = 0 ;

        N = 0;
        fid = fopen(fullfile(p,'doc','leadidtable_scpecg.txt'),'r');
        s = char(fread(fid,[1,inf],'uint8')); 
        fclose(fid);
        
        Code = repmat(NaN, 200, 1); Phi = Code; Theta = Code;
        while ~isempty(s),
        	[t,s] = strtok(s,[10,13]);
                if ~length(t)
                elseif ~strncmp(t,'#',1)    
                	ix3 = strfind(t,'MDC_ECG_LEAD_');
                	if isempty(ix3)
                		ix3 = length(t)+1;
                	end
                        [t1,t2] = strtok(t(1:ix3-1),[9,32]);
                        [t2,t3] = strtok(t2,[9,32]);
                        id = str2double(t2);
                        N  = N + 1;
                        Labels{N,1}	  = t1;
                        Code(N,1)         = id;
                        Description{N,1}  = deblank(t3);
                        %MDC_ECG_LEAD{N,1} = t(ix3+13:end)
                        MDC_ECG_LEAD{N,1} = t(ix3:end);
                end;
        end;
        N1 = N;

        % load table 
        fid = fopen(fullfile(p,'doc','elecpos.txt'),'r');
        t = char(fread(fid,[1,inf],'uint8'));
        fclose(fid);

        % extract table information       
        while ~isempty(t)
                [x,t] = strtok(t,[10,13]);
                if isempty(x)
                elseif strncmp(x,'#',1)
                else
                        N = N + 1;
                        [num,status,strarray] = str2double(x);
                        Code(N,1)   = num(1);
                        Labels{N,1} = upper(strarray{2});
                        Phi(N,1)    = num(3);
                        Theta(N,1)  = num(4);
                end;
        end;
        Phi   = Phi(:)  *pi/180;
        Theta = Theta(:)*pi/180;

        
        % loading is done only once. 
        BIOSIG_GLOBAL.XYZ = [sin(Theta).*cos(Phi), sin(Theta).*sin(Phi), cos(Theta)];
        BIOSIG_GLOBAL.Phi          = Phi*180/pi;
        BIOSIG_GLOBAL.Theta        = Theta*180/pi;
        BIOSIG_GLOBAL.LeadIdCode   = Code;
        BIOSIG_GLOBAL.Label        = Labels;
        BIOSIG_GLOBAL.Description  = Description;
        BIOSIG_GLOBAL.MDC_ECG_LEAD = MDC_ECG_LEAD;

        BIOSIG_GLOBAL.ISLOADED_XYZ = 1;
end; 


if nargin<1,
        HDR.LeadIdCode = BIOSIG_GLOBAL.LeadIdCode;
        HDR.Label      = BIOSIG_GLOBAL.Label;
        HDR.ELEC.XYZ   = BIOSIG_GLOBAL.XYZ; 
        HDR.TYPE       = 'ELPOS'; 
        
else    % electrode code and position

        if isstruct(arg1)
                HDR = arg1; 
        elseif isnumeric(arg1),
                HDR.LeadIdCode = arg1; 
        else
                HDR.Label = arg1; 
        end;

        tmp.flag1 = isfield(HDR,'ELEC');
        if tmp.flag1,
                tmp.flag1 = isfield(HDR.ELEC,'XYZ');
        end;
        if tmp.flag1,
                tmp.flag1 = any(HDR.ELEC.XYZ(:));
        end;
        tmp.flag2 = isfield(HDR,'LeadIdCode');
        tmp.flag3 = isfield(HDR,'Label');

        if (~tmp.flag1 || ~tmp.flag2 || ~tmp.flag3),
        	if 0, 
        	elseif tmp.flag3,
                        if ischar(HDR.Label)
                                HDR.Label = cellstr(HDR.Label);
                        end;
	                NS = length(HDR.Label); 
	        elseif tmp.flag2,
	        	NS = length(HDR.LeadIdCode); 
        	elseif isfield(HDR,'NS')
        		NS = HDR.NS;
        		HDR.LeadIdCode = zeros(1,HDR.NS);
	        end;      

                if tmp.flag3,
	                if ~tmp.flag1,
				HDR.ELEC.XYZ   = repmat(NaN,NS,3);
				HDR.ELEC.Phi   = repmat(NaN,NS,1);
				HDR.ELEC.Theta = repmat(NaN,NS,1);
	        	end;
                	if ~tmp.flag2,
                       		HDR.LeadIdCode = repmat(NaN,NS,1);
	        	end;
                	for k = 1:NS;
				Label = upper(deblank(HDR.Label{k})); 
				pos = find(Label==':');
				if ~isempty(pos)
					Label = Label(pos+1:end);
                            	end;        
	                        ix = strmatch(Label,BIOSIG_GLOBAL.Label,'exact');

	                        if length(ix)==2,
	                        	%%%%% THIS IS A HACK %%%%%
	                        	%% solve ambiguity for 'A1','A2'; could be EEG or ECG
	                        	if sum(HDR.LeadIdCode(1:k)>=996)>sum(HDR.LeadIdCode(1:k)<996)
	                        		%% majority are EEG electrodes,
	                        		ix = ix(find(BIOSIG_GLOBAL.LeadIdCode(ix)>996));
	                        	else	
	                        		%% majority are ECG electrodes,
	                        		ix = ix(find(BIOSIG_GLOBAL.LeadIdCode(ix)<996));
	                        	end;
	                        elseif isempty(ix)	
		                        ix = strmatch(deblank(HDR.Label{k}),BIOSIG_GLOBAL.MDC_ECG_LEAD,'exact');
	                        end; 	

        	                if (length(ix)==1),
	                	        if ~tmp.flag1,
        	                	        HDR.ELEC.XYZ(k,1:3) = BIOSIG_GLOBAL.XYZ(ix,:);
                	        	        HDR.ELEC.Phi(k)   = BIOSIG_GLOBAL.Phi(ix);
                        		        HDR.ELEC.Theta(k) = BIOSIG_GLOBAL.Theta(ix);
                        		end;
                	        	if ~tmp.flag2,
                                		HDR.LeadIdCode(k,1) = BIOSIG_GLOBAL.LeadIdCode(ix);
	                        	end;
        	                end;
                        end;
                else
			HDR.Label = cell(NS,1);
			for k = 1:NS;
				ix = find(BIOSIG_GLOBAL.LeadIdCode==HDR.LeadIdCode(k));
				if (length(ix)>=1),
					ix=ix(1);	% 
	                                HDR.Label{k} = BIOSIG_GLOBAL.Label{ix};
					if ~tmp.flag1,
        		                        HDR.ELEC.XYZ(k,1:3) = BIOSIG_GLOBAL.XYZ(ix,1:3);
        		                        HDR.ELEC.Phi(k,1)   = BIOSIG_GLOBAL.Phi(ix);
        		                        HDR.ELEC.Theta(k,1) = BIOSIG_GLOBAL.Theta(ix);
					end;
				else
					HDR.Label{k} = ['#',int2str(k)];
				end;
			end;
                end;
	end;
	if tmp.flag3 && ~any(HDR.LeadIdCode),
               	for k = 1:HDR.NS;
			ix = strmatch(upper(HDR.Label{k}),BIOSIG_GLOBAL.Label,'exact'); 
			if length(ix)==1,
				HDR.LeadIdCode(k) = BIOSIG_GLOBAL.LeadIdCode(ix);
			end;
		end; 		
        end;
end;

