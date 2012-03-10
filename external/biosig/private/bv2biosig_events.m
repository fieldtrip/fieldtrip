function HDR=bv2biosig_events(EVENT)
% BV2BIOSIG_EVENTS converts VMRK marker information BioSig Event codes. 
%  according to biosig/doc/eventcodes.txt. 
%  Currently, the convention of the BerlinBCI is implemented and supported. 
%
%  HDR = bv2biosig_events(arg1)   
%
%  arg1 can be an HDR-struct containg  HDR.EVENT.Desc 
%  or a struct containing EVENT.Desc 
%  or a cell-array Desc
%  or a char-array Desc 
% 
%  Warning: Approximately 32 (out of 510) events are currently not supported. 
%  For your data, you can check this with this command: 
%  HDR.EVENT.Desc(HDR.EVENT.TYP==0) 
% 
% see also: doc/eventcodes.txt

%	$Id$
%	Copyright (C) 2006,2007 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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


HDR.EVENT.Desc = []; 
if isfield(EVENT,'EVENT')
	HDR = EVENT; 
elseif isfield(EVENT,'Desc')
	HDR.EVENT = EVENT; 
elseif iscell(EVENT)
	HDR.EVENT.Desc = EVENT; 
elseif ischar(EVENT)
	for k=1:size(EVENT,1),
		HDR.EVENT.Desc{k} = EVENT(k,:); 
	end; 
else 
	error('unknown input argument')
end; 	

if isfield(HDR.EVENT,'TeegType')
	ix = strmatch('New Segment',HDR.EVENT.TeegType); 
	HDR.EVENT.TYP(ix)=hex2dec('7ffe'); 
end; 

if ~isfield(HDR,'NS')
	HDR.NS = NaN; 
end; 
if ~isfield(HDR,'Label'),
	HDR.Label = repmat({' '},min(HDR.NS,128),1);
end; 

FLAG_SEASON2_ARTERAWDATA = 0; 
if (HDR.NS==128)
	tmp = strvcat(HDR.Label);
	tmp(55:56,4) = 'f';
	FLAG_SEASON2_ARTERAWDATA = isequal(tmp(1:64,1:end-1),tmp(65:128,2:end)) & all(tmp(65:128,1)=='x');
	FLAG_SEASON2_ARTERAWDATA = FLAG_SEASON2_ARTERAWDATA & strncmp(HDR.FILE.Name,'arte',4);
end;
FLAG_ARTERAWDATA = strncmp(HDR.FILE.Name,'arte',4);

for k1 = 1:length(HDR.EVENT.POS)
	if isfield(HDR.EVENT,'Desc')
		tmp = HDR.EVENT.Desc{k1};
	elseif isfield(HDR.EVENT,'CodeDesc') && (HDR.EVENT.TYP(k1)<=length(HDR.EVENT.CodeDesc))
		tmp = HDR.EVENT.CodeDesc{HDR.EVENT.TYP(k1)};
	else 
		continue; 
	end;

	if 0,

	elseif strncmp(tmp,'TargetCode',10)
		HDR.EVENT.TYP(k1) = str2double(tmp(11:12))+hex2dec('0300'); 

	elseif strcmp(tmp,'BeginOfTrial')
		HDR.EVENT.TYP(k1) = hex2dec('0300'); 

	elseif strcmp(tmp,'hit')
        	HDR.EVENT.TYP(k1) = hex2dec('0381'); 
        elseif strcmp(tmp,'wrong')
        	HDR.EVENT.TYP(k1) = hex2dec('0382'); 

% eye movements
        elseif strcmpi(tmp,'augen links')
        	HDR.EVENT.TYP(k1) = hex2dec('0431');
        elseif strcmpi(tmp,'augen rechts')
        	HDR.EVENT.TYP(k1) = hex2dec('0432');
        elseif strcmpi(tmp,'augen hoch') | strcmpi(tmp,'augen oben')
        	HDR.EVENT.TYP(k1) = hex2dec('0433');
        elseif strcmpi(tmp,'augen unten') | strcmpi(tmp,'augen runter')
        	HDR.EVENT.TYP(k1) = hex2dec('0434');
        elseif strcmpi(tmp,'augen offen') | strcmp(tmp,'Augen offen & entspannen')
        	HDR.EVENT.TYP(k1) = hex2dec('0114');
        elseif strcmpi(tmp,'augen zu') | strcmp(tmp,'Augen zu & entspannen')
        	HDR.EVENT.TYP(k1) = hex2dec('0115');
        elseif strcmp(tmp,'blinzeln')
        	HDR.EVENT.TYP(k1) = hex2dec('0439'); 

% muscle movements 
        elseif strcmp(tmp,'EMG links')
        	HDR.EVENT.TYP(k1) = hex2dec('0441'); 
        elseif strcmp(tmp,'EMG rechts')
        	HDR.EVENT.TYP(k1) = hex2dec('0442'); 
        elseif strcmpi(tmp,'kopf bewegen')
        	HDR.EVENT.TYP(k1) = hex2dec('0443'); 
        elseif strcmp(tmp,'zunge an')
        	HDR.EVENT.TYP(k1) = hex2dec('0444'); 
        elseif strcmp(tmp,'zunge aus')
        	HDR.EVENT.TYP(k1) = hex2dec('8444'); 
        elseif strcmp(tmp,'Kiefer anspannen')
        	HDR.EVENT.TYP(k1) = hex2dec('0446'); 
        elseif strcmp(tmp,'beißen') | strcmp(tmp,'beissen') | strcmp(tmp,['bei',223,'en']),
        	HDR.EVENT.TYP(k1) = hex2dec('0446'); 
        elseif strcmp(tmp,'EMG fuss')
		HDR.EVENT.TYP(k1) = hex2dec('0447'); 
	elseif strcmp(tmp,'Arme bewegen')
		HDR.EVENT.TYP(k1) = hex2dec('0449'); 

% encoding der season2-arte* rawdata records
	elseif strncmp(tmp,'S',1) & (FLAG_SEASON2_ARTERAWDATA | FLAG_ARTERAWDATA), 
		n = str2double(tmp(2:end));
		if n==11,	% EMG left
			HDR.EVENT.TYP(k1) = hex2dec('0441'); 
		elseif n==12,	% EMG right
			HDR.EVENT.TYP(k1) = hex2dec('0442'); 
		elseif n==13,	% EMG (foot)
			HDR.EVENT.TYP(k1) = hex2dec('0447'); 
		elseif n==1,	% Augen (left)
			HDR.EVENT.TYP(k1) = hex2dec('0431'); 
		elseif n==2,	% Augen (right)
			HDR.EVENT.TYP(k1) = hex2dec('0432'); 
		elseif n==3,	% Augen oben
			HDR.EVENT.TYP(k1) = hex2dec('0433');
		elseif n==4,	% Augen unten
			HDR.EVENT.TYP(k1) = hex2dec('0434');
		elseif n==5,	% blinzeln
			HDR.EVENT.TYP(k1) = hex2dec('0439');
		elseif n==6,	% Augen zu & entspannen
			HDR.EVENT.TYP(k1) = hex2dec('0115');
		elseif n==7,	% Augen offen & entspannen
			HDR.EVENT.TYP(k1) = hex2dec('0114');
		elseif n==8,	% beissen
			HDR.EVENT.TYP(k1) = hex2dec('0446');
		elseif n==9,	% kopf bewegen
			HDR.EVENT.TYP(k1) = hex2dec('0443');
		elseif any(n==[10,100])	% ende der aktion
			HDR.EVENT.TYP(k1) = bitxor(hex2dec('8000'),HDR.EVENT.TYP(k1-1)); 
		else
			HDR.EVENT.TYP(k1) = n; 
		end;

% hits and misses, feedback
	elseif strncmp(tmp,'S',1) | strncmp(tmp,'R',1) 
		HDR.EVENT.CHN(k1) = (tmp(1)=='R')*64+1;  %%% hack to distinguish Player 1 and 2 in SEASON2 data
		n = str2double(tmp(2:end)); 
		if n==11,	% hit (left)
		       	HDR.EVENT.TYP(k1) = hex2dec('0381'); 
		elseif n==12,	% hit (right)
	        	HDR.EVENT.TYP(k1) = hex2dec('0381'); 
		elseif n==13,	% hit (foot)
	        	HDR.EVENT.TYP(k1) = hex2dec('0381'); 
        	elseif n==21,	% miss (left)
	        	HDR.EVENT.TYP(k1) = hex2dec('0382'); 
		elseif n==22,	% miss (right)
	        	HDR.EVENT.TYP(k1) = hex2dec('0382'); 
		elseif n==23,	% time out
	        	HDR.EVENT.TYP(k1) = hex2dec('830d');
		elseif n==60,	% feedback onset
	        	HDR.EVENT.TYP(k1) = hex2dec('030d'); 
		elseif any(n==[4,5,7]) %%% ignore these 
	        	HDR.EVENT.TYP(k1) = NaN; 
        	else
	        	HDR.EVENT.TYP(k1) = n; 
		end; 

% end of segment        	
        elseif strcmp(tmp,'s') | strcmp(tmp,'stop') | strcmp(tmp,'stopp'),
        	HDR.EVENT.TYP(k1) = bitxor(hex2dec('8000'),HDR.EVENT.TYP(k1-1)); 

        elseif ~isempty(tmp)
        	[n,v,s] = str2double(tmp(2:end));
        	if (length(n)==1) & (~v)
        		HDR.EVENT.TYP(k1) = n; 
       		end; 
        end; 	
end; 
HDR.EVENT.TYP = HDR.EVENT.TYP(:); 

if isfield(HDR.EVENT,'POS'); 
	% remove S4, S5, S7 before computing trial end
	ix = find(HDR.EVENT.TYP(2:end)==hex2dec('7ffe'))+1;
	for k = 1:length(ix), 
		if any(HDR.EVENT.TYP(ix(k)-1)==[1:3])
			HDR.EVENT.TYP(ix(k)-1) = NaN;	% if cue before end of segment is 1,2,or 3. 
		end; 	
	end; 
	
       	flag_remove = isnan(HDR.EVENT.TYP);
	HDR.EVENT.TYP = HDR.EVENT.TYP(~flag_remove);
	HDR.EVENT.POS = HDR.EVENT.POS(~flag_remove);
	HDR.EVENT.CHN = HDR.EVENT.CHN(~flag_remove);
	HDR.EVENT.DUR = HDR.EVENT.DUR(~flag_remove);
	if isfield(HDR.EVENT,'TeegType');
		HDR.EVENT.TeegType = HDR.EVENT.TeegType(~flag_remove);
	end;	
	if isfield(HDR.EVENT,'Desc')
		HDR.EVENT.Desc = HDR.EVENT.Desc(~flag_remove);
	end; 	

       	ix1 = find(HDR.EVENT.TYP<10);
       	ix2 = find(HDR.EVENT.TYP==100); 
       	if ~isempty(ix2),
		HDR.EVENT.TYP(ix2,1) = HDR.EVENT.TYP(ix2-1) + hex2dec('8000'); 
	end; 	
	ix0 = find((HDR.EVENT.TYP>0)&(HDR.EVENT.TYP<10));

	HDR.TRIG = HDR.EVENT.POS(ix0); 
	HDR.Classlabel = HDR.EVENT.TYP(ix0); 
end; 

if any(HDR.EVENT.DUR~=1)
	warning('Duration is not 1')
end; 

% convert from Type1 into Type3 table.
if 1, % ~isfield(HDR.EVENT,'CHN') & ~isfield(HDR.EVENT,'DUR'),  
	% HDR.EVENT.CHN = zeros(size(HDR.EVENT.POS)); 
	HDR.EVENT.DUR = zeros(size(HDR.EVENT.POS)); 

	% convert EVENT.Version 1 to 3, currently used by GDF, BDF and alpha
	flag_remove = zeros(size(HDR.EVENT.TYP));
	types  = unique(HDR.EVENT.TYP);
	for k1 = find(bitand(types(:)',hex2dec('8000')));
	        TYP0 = bitand(types(k1),hex2dec('7fff'));
	        TYP1 = types(k1);
	        ix0  = (HDR.EVENT.TYP==TYP0);
	        ix1  = (HDR.EVENT.TYP==TYP1);

	        if sum(ix0)==sum(ix1), 
	                HDR.EVENT.DUR(ix0) = HDR.EVENT.POS(ix1) - HDR.EVENT.POS(ix0);
	                flag_remove = flag_remove | (HDR.EVENT.TYP==TYP1);
                elseif 0, 
                	
                else 
	                fprintf(2,'Warning BV2BIOSIG_EVENT: number of event onset (TYP=%s) and event offset (TYP=%s) differ (%i-%i) in %s\n',dec2hex(double(TYP0)),dec2hex(double(TYP1)),sum(ix0),sum(ix1),HDR.FileName);
                        %% double(.) operator needed because Matlab6.5 can not fix(uint16(..))
	        end;
	end;
	if any(HDR.EVENT.DUR<0)
	        fprintf(2,'Warning SOPEN: EVENT ONSET later than EVENT OFFSET\n',dec2hex(TYP0),dec2hex(TYP1));
	        %HDR.EVENT.DUR(:) = 0
	end;
	HDR.EVENT.TYP = HDR.EVENT.TYP(~flag_remove);
	HDR.EVENT.POS = HDR.EVENT.POS(~flag_remove);
	HDR.EVENT.CHN = HDR.EVENT.CHN(~flag_remove);
	HDR.EVENT.DUR = HDR.EVENT.DUR(~flag_remove);
	if isfield(HDR.EVENT,'TeegType');
		HDR.EVENT.TeegType = HDR.EVENT.TeegType(~flag_remove);
	end;	
	%HDR.EVENT.TeegType = HDR.EVENT.TeegType(~flag_remove);
	if isfield(HDR.EVENT,'Desc')
		HDR.EVENT.Desc = HDR.EVENT.Desc(~flag_remove);
	end;
end;	
