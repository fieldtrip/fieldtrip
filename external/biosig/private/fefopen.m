function [HDR]=fefopen(arg1)
% FEFOPEN opens and reads FEF file 
%
%       HDR=fefopen(HDR)
%
%
% see also: SOPEN 
%
% References: 
% <A HREF="ftp://sigftp.cs.tut.fi/pub/eeg-data/standards/cenf060.zip ">About CEN/TC251</A> 

%	$Id$
%	Copyright (c) 1998,2004,2008  Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.

if ischar(arg1)
        HDR.FILE.OPEN = 0; 
        HDR.FileName = arg1;
        HDR.FILE.FID = fopen(HDR.FileName,'r');
        HDR.FILE.OPEN = 0; 
else
	HDR = arg1;
end;    


TYP_OID16 = 'uint16';

N = 0; 
tag=fread(HDR.FILE.FID,1,TYP_OID16);
len=fread(HDR.FILE.FID,1,'uint16');
while ~feof(HDR.FILE.FID),
        value=fread(HDR.FILE.FID,len,'uint8');
	
	tmp.TAG = sprintf('%x',tag); 
	tmp.LEN = len;
	tmp.VAL = value(1:min(len,10))';

        if tag>0, 
		N = N+1; 
		HDR.tag{N} = tmp;
	end;	 	
        
        tag=fread(HDR.FILE.FID,1,TYP_OID16);
        len=fread(HDR.FILE.FID,1,'uint16');
end;
        

if ~HDR.FILE.OPEN, 
        fclose(HDR.FILE.FID);
end; 
