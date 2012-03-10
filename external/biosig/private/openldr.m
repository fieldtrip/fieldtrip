function LDR = openldr(FN,PERMISSION,Mode,arg4,arg5,arg6)
% OPENLDR loads neuroscan LDR files
% LDR = OPENLDR(Filename [, PERMISSION [, Mode]]);
% LDR = OPENLDR(LDR [, PERMISSION [, Mode]]);
%
% LDR is a struct with the following fields
%   LDR.FileName 	Name of LDR-file
%   LDR.Label_Out	Labels of output channels
%   LDR.Label_In	Labels of input channels
%   LDR.RR		re-referencing matrix
%   LDR.datatype	'REREF_MATRIX' indicates this datatype
%
% PERMISSION	'r'	reads LDR file
% 		'w'	writes LDR file
% 		'r+w'	reads and writes LDR file 
%			(useful in combination with RESCALE-Mode) 
%
% Mode [optional] 'RESCALE' performs a rescaling of the weights
%		sum of positive weights becomes +1
%		sum of negative weights becomes -1
%

%	$Id$
%	Copyright (C) 1997-2003,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.


if nargin<2, PERMISSION=''; end;
if nargin<3, Mode=''; end;

if ~isstruct(FN), 
        LDR.FileName = FN;
        PERMISSION = [PERMISSION,'r']; 
else
        LDR = FN;
        PERMISSION = [PERMISSION,'w']; 
end;


if any(PERMISSION=='r'); 
        % load LDR-file 
        fid = fopen(LDR.FileName);
        if fid<0, fprintf(2,'File %s not found.\n',LDR.FileName); return; end; 
        
        % reads line 1: size information
        s  = fgetl(fid);
        if ~ischar(s), fprintf(2,'ERROR LOADLDR: file %s corrupted\n',FN); end;
        sz = str2num(s);
        
        % reads line 2: output labels
        s  = fgetl(fid);
        if ~ischar(s), fprintf(2,'ERROR LOADLDR: file %s corrupted\n',FN); end;
        tmp = reshape(s,12,sz(2)+1)';
        LDR.Label_Out = tmp(2:sz(2)+1,:);
        
        % read lines 3+: input labels and weights
        for k = 1:sz(1),
                s = fgetl(fid);
                if ~ischar(s), fprintf(2,'ERROR LOADLDR: file %s corrupted\n',FN); end;
                r = reshape(s,12,sz(2)+1)';
                LDR.Label_In(k,1:12) = char(abs(s(1:12)));
                LDR.RR(k,1:sz(2)) = str2num(r(2:size(r,1),:))';
        end;
        
        fclose(fid);
        LDR.datatype='REREF_MATRIX';
        
end;


if strcmp(Mode,'RESCALE');
        tmp = LDR.RR;
        tmp(tmp<0) = 0;
        w = sum(tmp);
        RR = zeros(sz);
        if ~all(w==0 | w==1)
                ix = find(w);
                RR(:,ix) = tmp(:,ix)./w(ones(sz(1),1),ix);
                % RR(:,ix) = (tmp(:,ix)>0)./(ones(sz(1),1)*sum(tmp(:,ix)>0));
        end;
        tmp = LDR.RR;
        tmp(tmp>0)=0;
        w = sum(tmp);
        if ~all(w==0 | w==-1)
                ix = find(w);
                RR(:,ix) = tmp(:,ix)./w(ones(sz(1),1),ix);
                % RR(:,ix) = RR(:,ix)-(tmp(:,ix)<0)./(ones(sz(1),1)*sum(tmp(:,ix)<0));
        end;
        LDR.RR=RR;
end;

if any(PERMISSION=='w'); 
        tmp = isfield(LDR,'datatype') & strcmp(LDR.datatype,'REREF_MATRIX');
        if ~tmp,
                warning('LDR.datatype does not fit');
        end;
        sz = size(LDR.RR);
        if sz(1)~=size(LDR.Label_In,1);
                fprintf(2,'Size of Label_In does not fit RR\n');
                return;
        end;
        if sz(2)~=size(LDR.Label_Out,1);
                fprintf(2,'Size of Label_Out does not fit RR\n');
                return;
        end;
        
	% open LDR-file 
        fid = fopen(LDR.FileName,'w+');
	if fid<0, fprintf(2,'Couldnot open file %s .\n',LDR.FileName); return; end; 
        
        % write line 1: size information
        fprintf(fid,'%i %i\n',sz);
        
        % condition label information
        %LDR.Label_Out = char(LDR.Label_Out);

        nc = size(LDR.Label_Out,2);
        if nc<12,
                LDR.Label_Out = [LDR.Label_Out,abs(' ')*ones(sz(1),12-nc)];
        elseif nc>12,
                LDR.Label_Out = LDR.Label_Out(:,1:12);
        end;
        %LDR.Label_In = char(LDR.Label_In);
        nc = size(LDR.Label_In,2);
        if nc<12,
                LDR.Label_In = [LDR.Label_In,abs(' ')*ones(sz(1),12-nc)];
        elseif nc>12,
                LDR.Label_In = LDR.Label_In(:,1:12);
        end;
        
	% write line 2: out-labels
        fwrite(fid,32+zeros(12,1),'uint8');
        %fprintf(fid,'%c',abs(' ')*ones(12,1));
	for k = 1:sz(2),
                fwrite(fid,abs(LDR.Label_Out(k,1:12)),'uint8');
        end;
        
        % write lines 3+: in-labels and weights
        for k = 1:sz(1),
                fprintf(fid,'\n');
                fwrite(fid,abs(LDR.Label_In(k,1:12)),'uint8');
                fprintf(fid,'%12.5f',LDR.RR(k,:));
        end;
        fclose(fid);
end;
