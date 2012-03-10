function [data,HDR] = iread(HDR,CHAN,StartPos)
% IREAD reads image data 
%
% [S,HDR] = iread(HDR)
%
% See also: fread, SREAD, SWRITE, SCLOSE, SSEEK, SREWIND, STELL, SEOF

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.

%	$Id$
%	(C) 2005,2006,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


if 0, 


elseif strcmp(HDR.TYPE,'IMAGE:BMP'),
        fseek(HDR.FILE.FID,HDR.HeadLen,'bof');
        nc = ceil((HDR.bits*HDR.IMAGE.Size(1))/32)*4;
        
        if (HDR.bits==1)
                [data,c] = fread(HDR.FILE.FID,[nc,HDR.IMAGE.Size(2)*8],'ubit1');
                %[data,c] = fread(HDR.FILE.FID,HDR.IMAGE.Size([2,1]),'ubit1');
size(data),
%                data = reshape(data(1:nc,:)',HDR.IMAGE.Size);
                
        elseif (HDR.bits==4)
                palr   = [  0,128,  0,128,  0,128,  0,192,128,255,  0,255,  0,255,  0,255]; 
                palg   = [  0,  0,128,128,  0,  0,128,192,128,  0,255,255,  0,  0,255,255]; 
                palb   = [  0,  0,  0,  0,128,128,128,192,128,  0,  0,  0,255,255,255,255]; 
                tmp    = uint8(fread(HDR.FILE.FID,[nc,HDR.IMAGE.Size(2)*2],'ubit4'));
                data        = palr(tmp(1:HDR.IMAGE.Size(1),:)'+1);
                data(:,:,2) = palg(tmp(1:HDR.IMAGE.Size(1),:)'+1);
                data(:,:,3) = palb(tmp(1:HDR.IMAGE.Size(1),:)'+1);
                data = data(HDR.IMAGE.Size(2):-1:1,:,:);
                
        elseif (HDR.bits==8)
                pal = uint8(colormap*256);
                tmp = fread(HDR.FILE.FID,[nc,HDR.IMAGE.Size(2)],'uint8');
                data        = pal(tmp(1:HDR.IMAGE.Size(1),:)'+1,1);
                data(:,:,2) = pal(tmp(1:HDR.IMAGE.Size(1),:)'+1,2);
                data(:,:,3) = pal(tmp(1:HDR.IMAGE.Size(1),:)'+1,3);
                data = reshape(data,[HDR.IMAGE.Size([2,1]),3]);
                data = data(HDR.IMAGE.Size(2):-1:1,:,:);
                
        elseif (HDR.bits==24)
                data   = uint8(fread(HDR.FILE.FID,[nc,HDR.IMAGE.Size(2)],'uint8'));
                HDR.BMP.Red   = data((1:HDR.IMAGE.Size(1))*3,:)';
                HDR.BMP.Green = data((1:HDR.IMAGE.Size(1))*3-1,:)';
                HDR.BMP.Blue  = data((1:HDR.IMAGE.Size(1))*3-2,:)';
                data = HDR.BMP.Red;
                data(:,:,2) = HDR.BMP.Green;
                data(:,:,3) = HDR.BMP.Blue;
                data = data(HDR.IMAGE.Size(2):-1:1,:,:);
        else
                
        end;
    

elseif strcmp(HDR.TYPE,'IMAGE:FITS'),
	if nargin>1,
		KK = CHAN;
	else	
		[tmp, KK] = max(HDR.IMAGE_Size);   % select block
	end;
	
	for KK = 1:length(HDR.FITS),

	status = fseek(HDR.FILE.FID,HDR.HeadLen(KK),'bof');

	HDR.AS.bps = abs(HDR.FITS{KK}.BITPIX)/8;
	if HDR.FITS{KK}.BITPIX == 8,
		HDR.GDFTYP = 'uint8';
	elseif HDR.FITS{KK}.BITPIX == 16,
		HDR.GDFTYP = 'int16';
	elseif HDR.FITS{KK}.BITPIX == 32,
		HDR.GDFTYP = 'int32';
	elseif HDR.FITS{KK}.BITPIX == 64,
		HDR.GDFTYP = 'int64';
	elseif HDR.FITS{KK}.BITPIX == -32,
		HDR.GDFTYP = 'float32';
	elseif HDR.FITS{KK}.BITPIX == -64,
		HDR.GDFTYP = 'float64';
	else
		warning('IREAD (FITS{KK})');
	end;	
	
	if isfield(HDR.FITS{KK},'BZERO')	HDR.Off = HDR.FITS{KK}.BZERO;
	else					HDR.Off = 0;			end;		
	if isfield(HDR.FITS{KK},'BSCALE')	HDR.Cal = HDR.FITS{KK}.BSCALE;
	else					HDR.Cal = 1;			end;		
	if isfield(HDR.FITS{KK},'BUNIT'),	HDR.PhysDim = HDR.FITS{KK}.BUNIT;
	else					HDR.PhysDim = '[1]';		end;		
	if isfield(HDR.FITS{KK},'DATAMAX'),	HDR.PhysMax = HDR.FITS{KK}.DATAMAX;
	else					HDR.PhysMax = NaN;		end;
	if isfield(HDR.FITS{KK},'DATAMIN'),	HDR.PhysMin = HDR.FITS{KK}.DATAMIN;
	else					HDR.PhysMin = NaN;	 	end;

	if ~isfield(HDR.FITS{KK},'XTENSION')
		[data,c] = fread(HDR.FILE.FID, prod(HDR.IMAGE(KK).Size), HDR.GDFTYP);
		data = reshape(data,HDR.IMAGE(KK).Size) * HDR.Cal + HDR.Off;

	elseif strncmp(HDR.FITS{KK}.XTENSION,'TABLE',5)
		for k = 1:HDR.FITS{KK}.TFIELDS,
			f = ['TTYPE',int2str(k)];
			if isfield(HDR.FITS{KK},f)
				HDR.TABLE{KK}.Label(k,:) = getfield(HDR.FITS{KK},f);
			end;	
			tmp = getfield(HDR.FITS{KK},['TFORM',int2str(k)]);
			typ(k) = tmp(1);
			ix(k)  = getfield(HDR.FITS{KK},['TBCOL',int2str(k)]);
		end; 
		ix(k+1) = HDR.FITS{KK}.NAXIS1+1;

		sa = {};
		[data,c] = fread(HDR.FILE.FID, HDR.IMAGE(KK).Size, HDR.GDFTYP);
		data = data';
		s = repmat(NaN, HDR.FITS{KK}.NAXIS2, HDR.FITS{KK}.TFIELDS);
		for k = 1:HDR.FITS{KK}.TFIELDS,
			[s(:,k), status,sa(:,k)] = str2double(char(data(:,ix(k):ix(k+1)-1)));
		end;
		HDR.TABLE{KK} = s;
        
	elseif strncmp(HDR.FITS{KK}.XTENSION,'BINTABLE',8)
		sz = ones(1,HDR.FITS{KK}.TFIELDS);
		for k = 1:HDR.FITS{KK}.TFIELDS,
			f = ['TTYPE',int2str(k)];
			if isfield(HDR.FITS{KK},f)
				HDR.TABLE{KK}.Label(k,:) = getfield(HDR.FITS{KK},f);
			end;
			tmp = getfield(HDR.FITS{KK},['TFORM',int2str(k)]);
			ix  = min(find(tmp>'9'));
			HDR.FITS{KK}.TYP(k) = tmp(ix);
			tmp1 = str2double(tmp(1:ix-1));
			if ~isempty(tmp1), sz(k)=tmp1; end;

			if 0, 
			elseif tmp(ix)=='L', 	GDFTYP{k} = 'uint8';	% char T, F
			elseif tmp(ix)=='X', 	GDFTYP{k} = 'ubit1';	%ubit1
			elseif tmp(ix)=='B', 	GDFTYP{k} = 'uint8';	% uint8
			elseif tmp(ix)=='I', 	GDFTYP{k} = 'int16';	% int16
			elseif tmp(ix)=='J', 	GDFTYP{k} = 'int32';	% int32
			elseif tmp(ix)=='A', 	GDFTYP{k} = 'uchar';	% char
			elseif tmp(ix)=='E', 	GDFTYP{k} = 'float32';	% float32
			elseif tmp(ix)=='D', 	GDFTYP{k} = 'float64';	% double
			elseif tmp(ix)=='C', sz(k) = sz(k)*2;	GDFTYP{k} = 'float32';	% [1+i]*float32
			elseif tmp(ix)=='M', sz(k) = sz(k)*2;	GDFTYP{k} = 'float64';	% [1+i]*double
			elseif tmp(ix)=='P', 	GDFTYP{k} = 'uint64';	% uint64, array desc ? 
			else
			end;
		end; 

		for k2 = 1:HDR.FITS{KK}.NAXIS2,
		for k1 = 1:HDR.FITS{KK}.TFIELDS,
			[s,c] = fread(HDR.FILE.FID, [1,sz(k1)], GDFTYP{k1});

			if 0, 
			elseif HDR.FITS{KK}.TYP(k1)=='L',
				sig{k1}(k2,:) = (s=='T');
			elseif HDR.FITS{KK}.TYP(k1)=='A',
				sig{k1}(k2,:) = char(s);
			elseif any(HDR.FITS{KK}.TYP(k1)=='CM'),
				sig{k1}(k2,:) = [1,i]*reshape(s,2,sz(k1)/2);
			else
				sig{k1}(k2,:) = s;
			end;	
 		end;
		end;
		HDR.TABLE{KK} = sig;
        
	elseif strncmp(HDR.FITS{KK}.XTENSION,'IMAGE',5)
		[data,c] = fread(HDR.FILE.FID, prod(HDR.IMAGE(KK).Size), HDR.GDFTYP);
		
		[c,size(data),HDR.IMAGE(KK).Size]
		data = reshape(data,HDR.IMAGE(KK).Size) * HDR.Cal + HDR.Off;

	else
		[data,c] = fread(HDR.FILE.FID, prod(HDR.IMAGE(KK).Size), HDR.GDFTYP);
		data = reshape(data,HDR.IMAGE(KK).Size) * HDR.Cal + HDR.Off;

	end;
	end;


elseif strcmp(HDR.TYPE,'IMAGE:HGT'),
	data = fread(HDR.FILE.FID,HDR.IMAGE.Size,'int16')'; 


elseif strcmp(HDR.TYPE,'IMAGE:TIFF'),

	if ~isfield(HDR.TIFF,'StripOffset')
        if 1, 
                data = []; count = 0; 
                for k = 1:length(HDR.TIFF.TileOffset)
                        status = fseek(HDR.FILE.FID,HDR.TIFF.TileOffset(k),'bof');
                        [d,c]  = fread(HDR.FILE.FID,HDR.TIFF.TileByteCount(k),HDR.GDFTYP);
                        data = [data; d]; count = count + c;
                end;	
	else

		fprintf(HDR.FILE.stderr,'Warning IREAD (TIFF): TileOffset not supported yet.\n');
                data = repmat(NaN,HDR.IMAGE.Size); 
                count = 0; 
                nr = HDR.IMAGE.Size(1)./HDR.TIFF.TileLength;
                nc = HDR.IMAGE.Size(2)./HDR.TIFF.TileWidth;
                for k1 = 1:nr
                for k2 = 1:nc,
                        status = fseek(HDR.FILE.FID,HDR.TIFF.TileOffset(k1*nc-nc+k2),'bof');
                        [d,c]  = fread(HDR.FILE.FID,HDR.TIFF.TileByteCount(k1*nc-nc+k2),HDR.GDFTYP);
                        d = permute(reshape(d,[HDR.TIFF.SamplesPerPixel,HDR.TIFF.TileWidth,HDR.TIFF.TileLength]),[3,2,1]);
                        [size(d),size(data)],
                        data ((k1-1)*HDR.TIFF.TileLength+1:k1*HDR.TIFF.TileLength,(k2-1)*HDR.TIFF.TileWidth+1:k2*HDR.TIFF.TileWidth,:) = d; 
                        count = count + c;
                end;
                end;	
        end;

       
        elseif ~any(HDR.TIFF.StripOffset(1:end-1)-HDR.TIFF.StripOffset(2:end)+HDR.TIFF.StripByteCounts(1:end-1))
                status = fseek(HDR.FILE.FID,HDR.TIFF.StripOffset(1),'bof');
                [data,count] = fread(HDR.FILE.FID,sum(HDR.TIFF.StripByteCounts)*length(HDR.Bits),HDR.GDFTYP);
        else
                data = []; count = 0; 
                for k = 1:length(HDR.TIFF.StripOffset)
                        status = fseek(HDR.FILE.FID,HDR.TIFF.StripOffset(k),'bof');
                        [d,c]  = fread(HDR.FILE.FID,HDR.TIFF.StripByteCounts(k),HDR.GDFTYP);
                        data = [data; d]; count = count + c;
                end;	
        end;

        
        if HDR.TIFF.Compression==1, 
                
        elseif HDR.TIFF.Compression==32773,	% Unpack PackBits
                D = data;
                N = prod(HDR.IMAGE.Size);
                data = zeros(N,1);
                k1 = 1; k2 = 0; 
                while k2 < N; %length(D),
                        if (D(k1)<128),
                                data(k2+(1:D(k1)+1)) = D(k1+(1:D(k1)+1));
                                k2 = k2 + D(k1) + 1;
                                k1 = k1 + D(k1) + 2;
                        elseif (D(k1)==128),
                                k1 = k1 + 1;
                        elseif (D(k1)>128),
                                n  = 1 -  (D(k1) - 256);
                                data(k2+(1:n)) = D(k1+1);
                                k1 = k1 + 2;
                                k2 = k2 + n;
                        else
                                error('TIFF decompression error');
                        end;	
                end;
                
                if (k2~=length(data)) | (k1-1-length(D))
                        [k2,length(data),k1,length(D)],
                        warning('TIFF decompression not completed');
                        data = data(1:N);
                end;
        else
                fprintf(HDR.FILE.stderr,'Error IREAD (TIFF): decompression Mode=%i not supported.\n',HDR.TIFF.Compression);
        end;

        if prod(HDR.IMAGE.Size)<=length(data);
        	tmp = length(data)-prod(HDR.IMAGE.Size);
	        if tmp,
			fprintf(HDR.FILE.stderr,'Warning IREAD (TIFF): %i bytes to much.\n',tmp);		
        	end;
        	       
                data = permute(reshape(data(1:prod(HDR.IMAGE.Size)),HDR.IMAGE.Size([3,2,1])),[3,2,1]);
	        if HDR.IMAGE.Size(3)==1,
			data= squeeze(data);
	        end;
        else
                fprintf(HDR.FILE.stderr,'Warning IREAD (TIFF): size of data does not fit.\n');	
        end;	
        
        if isfield(HDR.TIFF,'ColorMap');
                HDR.IMAGE.ColorMap = HDR.TIFF.ColorMap*2^-16;
                
        elseif isfield(HDR.TIFF,'SamplesPerPixel');
                if HDR.TIFF.SamplesPerPixel==1,	
                        data = data*2^(-HDR.Bits);
                end;
        end;
    

end;
	
	
