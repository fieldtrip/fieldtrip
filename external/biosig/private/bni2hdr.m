function [HDR]=bni2hdr(arg1,arg3,arg4,arg5,arg6)
% BNI2HDR converts BNI header information into BioSig Header information
%	HDR = BNI2HDR(HDR)
%
% INPUT:
%   HDR.H1 contains ascii header 
%
% OUTPUT:
%   HDR.Label
%   HDR.T0 
%       ...
%
% see also: SOPEN 

%	$Id$
%	Copyright (c) 2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
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
elseif exist(arg1,'file');
	HDR.FileName = arg1; 
	[HDR.FILE.Path,HDR.FILE.Name,HDR.FILE.Ext]=fileparts(arg1);
end;
if ~isfield(HDR,'H1')
	fid = fopen(fullfile(HDR.FILE.Path,HDR.FILE.Name,'.bni'),'rt');
	HDR.H1 = char(fread(fid,[1,inf],'uint8')); 
	fclose(fid);
end;

s = HDR.H1;
if ~strcmp(s(1:28),'FileFormat = BNI-1-BALTIMORE')
	fprintf(HDR.FILE.stderr,'WARNING BNI2HDR: Header information is not Nicolet BNI format.\n'); 
end;	

while ~isempty(s),
	[t,s]=strtok(s,[10,13]);
	[t1,t2]=strtok(t,' =');  
	[t2,t3]=strtok(t2,' =');  
	if strcmp(t1,'PatientId')
		HDR.Patient.Id = t2;
	elseif strcmpi(t1,'Sex')
		HDR.Patient.Sex = t2; % strncmpi(t2,'m',1)+2*strncmpi(t2,'f',1);
	elseif strncmpi(t1,'medication',10)
		HDR.Patient.Medication = t2;
	elseif strncmpi(t1,'diagnosis',10)
		HDR.Patient.Diagnosis = t2;
	elseif strcmpi(t1,'MontageRaw')
		[tmp1,tmp2,HDR.Label] = str2double(t2,',');
	elseif strcmpi(t1,'Age')
		HDR.Patient.Age = str2double(t2);
	elseif strcmp(t1,'Date')
		if any(t2=='/')
			t2(t2=='/')=' ';
			HDR.T0([2,3,1])=str2double(t2); 	
		end; 
	elseif strcmp(t1,'Time')
		t2(t2==':') = ' ';
		HDR.T0(4:6) = str2double(t2);
	elseif strcmp(t1,'Rate')
		HDR.SampleRate = str2double(t2);
	elseif strcmp(t1,'NchanFile')
		HDR.NS = str2double(t2);
	elseif strcmp(t1,'UvPerBit')
		HDR.Cal = str2double(t2);
	elseif strcmp(t1,'[Events]')
		% not supported yet 
		s = [];
	end; 
end; 
