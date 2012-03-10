function [datatyp,limits,datatypes,numbits,GDFTYP] = gdfdatatype(GDFTYP)
% GDFDATATYPE converts number into data type according to the definition of the GDF format [1]. 
%
%   [datatyp,limits,datatypes,numbits,GDFTYP]=gdfdatatype(gdftyp)
%
%
% See also: SOPEN, SREAD, SWRITE, SCLOSE
%
% References:
% [1] A. Schlögl, O. Filz, H. Ramoser, G. Pfurtscheller, GDF - A general dataformat for biosignals, Technical Report, 2004.
% available online at: http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/eeg/gdf4/TR_GDF.pdf. 

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.

%	$Id$
%	(C) 1997-2005,2008 by Alois Schloegl <a.schloegl@ieee.org>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


if ischar(GDFTYP),
        GDFTYP = cellstr(GDFTYP); 
end;
if iscell(GDFTYP),
        n = length(GDFTYP); 
        gdftyp = repmat(NaN,size(GDFTYP)); 
        for k = 1:n; 
                if strncmpi(GDFTYP{k},'char',4)
                        gdftyp(k) = 0;
                elseif strncmpi(GDFTYP{k},'int8',4)
                        gdftyp(k) = 1;
                elseif strncmpi(GDFTYP{k},'uint8',5)
                        gdftyp(k) = 2;
                elseif strncmpi(GDFTYP{k},'int16',5)
                        gdftyp(k) = 3;
                elseif strncmpi(GDFTYP{k},'uint16',6)
                        gdftyp(k) = 4;
                elseif strncmpi(GDFTYP{k},'int32',5)
                        gdftyp(k) = 5;
                elseif strncmpi(GDFTYP{k},'uint32',6)
                        gdftyp(k) = 6;
                elseif strncmpi(GDFTYP{k},'int64',5)
                        gdftyp(k) = 7;
                elseif strncmpi(GDFTYP{k},'uint64',6)
                        gdftyp(k) = 8;
                elseif strncmpi(GDFTYP{k},'float32',7) || strcmp(GDFTYP{k},'float')
                        gdftyp(k) = 16;
                elseif strncmpi(GDFTYP{k},'float64',7) || strcmp(GDFTYP{k},'double')
                        gdftyp(k) = 17;
                elseif strncmpi(GDFTYP{k},'float128',7)
                        gdftyp(k) = 18;
                elseif strncmpi(GDFTYP{k},'bit',3)
                        [num,status] = str2double(GDFTYP{k}(4:end));
                        if ~status
                                gdftyp(k) = 255+num;
                        end;
                elseif strncmpi(GDFTYP{k},'ubit',4)
                        [num,status] = str2double(GDFTYP{k}(5:end));
                        if ~status
                                gdftyp(k) = 511+num;
                        end;
                end
        end;
        GDFTYP = gdftyp; 
end;

datatyp = [];
datatypes = {};
limits = repmat(NaN,length(GDFTYP),3);
numbits = repmat(NaN,size(GDFTYP));

for k=1:length(GDFTYP),
        if GDFTYP(k)==0
                datatyp=('uint8');
                limit = [0,256,256];       
                nbits = 8; 
        elseif GDFTYP(k)==1
                datatyp=('int8');
                limit = [-128,127,128];        
                nbits = 8; 
        elseif GDFTYP(k)==2
                datatyp=('uint8');
                limit = [0,256,256];        
                nbits = 8; 
        elseif GDFTYP(k)==3
                datatyp=('int16');
                limit = [-2^15,2^15-1,-2^15];        
                nbits = 16; 
        elseif GDFTYP(k)==4
                datatyp=('uint16');
                limit = [0,2^16,2^16];        
                nbits = 16; 
        elseif GDFTYP(k)==5
                datatyp=('int32');
                limit = [-2^31,2^31-1,-2^31];        
                nbits = 32; 
        elseif GDFTYP(k)==6
                datatyp=('uint32');
                limit = [0,2^32,2^32];        
                nbits = 32; 
        elseif GDFTYP(k)==7
                datatyp=('int64');
                limit = [-2^63,2^63-1,-2^63];        
                nbits = 64; 
        elseif GDFTYP(k)==8
                datatyp=('uint64');
                limit = [0,2^64,2^64];        
                nbits = 64; 
        elseif GDFTYP(k)==16
                datatyp=('float32');
                limit = [-(2-2^-23)*2^127,(2-2^-23)*2^127,NaN];        
                nbits = 32; 
        elseif GDFTYP(k)==17
                datatyp=('float64');
                limit = [-realmax,realmax,NaN];        
                nbits = 64; 
        elseif GDFTYP(k)==18
                datatyp=('float128');
                limit = [-inf,inf,NaN];        
                nbits = 128; 
        elseif (GDFTYP(k)>255) && (GDFTYP(k)<512)
                nbits = GDFTYP(k)-255;
                datatyp=['bit',int2str(nbits)];
                limit = [-(2^(nbits-1)),2^(nbits-1)-1,-(2^(nbits-1))];
        elseif (GDFTYP(k)>511) && (GDFTYP(k)<768)
                nbits = GDFTYP(k)-511;
                datatyp=['ubit',int2str(nbits)];
                limit = [0,2^nbits,2^nbits];
        else 
                fprintf(2,'Error GDFREAD: Invalid GDF channel type\n');
                datatyp='';
                limit = [NaN,NaN,NaN];
                nbits = NaN; 
        end;
        datatypes{k}  = datatyp;
        limits(k,:)   = limit;
        numbits(k)    = nbits; 
end;