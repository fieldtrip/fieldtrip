function T = spm_type(x, arg)
% translates data type specifiers between SPM & Matlab representations
% FORMAT T = spm_type(x, arg)
% x    - specifier
% T    - type
% arg  - optional string argument, can be
%	 - 'swapped' - if type is byteswapped return 1.
%	 - 'maxval'  - return maximum allowed value.
%	 - 'minval'  - return minimum allowed value.
%	 - 'nanrep'  - return 1 if there is a NaN representation.
%	 - 'bits'    - return the number of bits per voxel.
%	 - 'intt'    - return 1 if values rounded to nearest integer.
%_______________________________________________________________________
%
% Original format specifiers are based on ANALYZE.  If the input is
% a number then the corresponding matlab string is returned by default.
% If the input is a string then the appropriate TYPE is returned.
% However, if the optional arg argument is supplied then other
% information will be returned instead.
%
% With no arguments, a list of data types is returned.
%
% Additional support was added for signed bytes, unsigned short and
% unsigned int (by adding 128 to the format specifiers for unsigned bytes
% signed short and signed int).  Byte swapped datatypes have the same
% identifiers as the non-byte-swapped versions, multiplied by a factor of
% 256.
%_______________________________________________________________________
% @(#)spm_type.m	2.3 John Ashburner, Andrew Holmes 99/04/27


prec = str2mat('uint8','int16','int32','float','double','int8','uint16','uint32','uint8','int16','int32','float','double','int8','uint16','uint32');
types   = [    2      4      8   16   64   130    132    136,   512   1024   2048 4096 16384 33280  33792  34816];
swapped = [    0      0      0    0    0     0      0      0,     1      1      1    1     1     1      1      1];
maxval  = [2^8-1 2^15-1 2^31-1  Inf  Inf 2^7-1 2^16-1 2^32-1, 2^8-1 2^15-1 2^31-1  Inf   Inf 2^8-1 2^16-1 2^32-1];
minval  = [    0  -2^15  -2^31 -Inf -Inf  -2^7      0      0,     0  -2^15  -2^31 -Inf  -Inf  -2^7      0      0];
nanrep  = [    0      0      0    1    1     0      0      0,     0      0      0    1     1     0      0      0];
bits    = [    8     16     32   32   64     8     16     32,     8     16     32   32    64     8     16     32];
intt    = [    1      1      1    0    0     1      1      1,     1      1      1    0     0     1      1      1];

if nargin==0,
	T=types;
	return;
end;

if ischar(x),
	sel = [];
	msk = find(swapped==0);
	for i=msk,
		if strcmp(deblank(prec(i,:)),deblank(x)), 
			sel = i;
			break;
		end;
	end;
else,
	sel = find(types == x);
end;
if nargin == 1,
	if ischar(x),
		if isempty(sel), T = NaN;
		else, T = types(sel); end;
	else,
		if isempty(sel), T = 'unknown';
		else, T = deblank(prec(sel,:)); end;
	end;
elseif isempty(sel),
	T = NaN;
else,
	switch lower(arg)
	case 'swapped', T = swapped(sel);
	case 'maxval',  T = maxval(sel);
	case 'minval',  T = minval(sel);
	case 'nanrep',  T = nanrep(sel);
	case 'bits',    T = bits(sel);
	case 'intt',    T = intt(sel);
	otherwise,      T = NaN;
	end;
end;
