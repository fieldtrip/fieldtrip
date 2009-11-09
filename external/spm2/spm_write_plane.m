function V = spm_write_plane(V,A,p)
% Write a transverse plane of image data.
% FORMAT V = spm_write_plane(V,A,p)
% V   - data structure containing image information.
%       - see spm_vol for a description.
% A   - the two dimensional image to write.
% p   - the plane number (beginning from 1).
%
% VO  - (possibly) modified data structure containing image information.
%       It is possible that future versions of spm_write_plane may
%       modify scalefactors (for example).
%
%_______________________________________________________________________
% @(#)spm_write_plane.m	2.19 John Ashburner 03/07/16

if any(V.dim(1:2) ~= size(A)), error('Incompatible image dimensions');      end;
if p>V.dim(3),                 error('Plane number too high');              end;

% Write Analyze image by default
V = write_analyze_plane(V,A,p);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function V = write_analyze_plane(V,A,p)

types   = [    2      4      8   16   64   130    132    136,   512   1024   2048 4096 16384 33280  33792  34816];
maxval  = [2^8-1 2^15-1 2^31-1  Inf  Inf 2^7-1 2^16-1 2^32-1, 2^8-1 2^15-1 2^31-1  Inf   Inf 2^8-1 2^16-1 2^32-1];
minval  = [    0  -2^15  -2^31 -Inf -Inf  -2^7      0      0,     0  -2^15  -2^31 -Inf  -Inf  -2^7      0      0];
intt    = [    1      1      1    0    0     1      1      1,     1      1      1    0     0     1      1      1];
prec = str2mat('uint8','int16','int32','float','double','int8','uint16','uint32','uint8','int16','int32','float','double','int8','uint16','uint32');
swapped = [    0      0      0    0    0     0      0      0,     1      1      1    1     1     1      1      1];
bits    = [    8     16     32   32   64     8     16     32,     8     16     32   32    64     8     16     32];

dt      = find(types==V.dim(4));
if isempty(dt), error('Unknown datatype'); end;

A = double(A);

% Rescale to fit appropriate range
if intt(dt),
	A(isnan(A)) = 0;
	mxv         = maxval(dt);
	mnv         = minval(dt);
	A           = round(A*(1/V.pinfo(1)) - V.pinfo(2));
	A(A > mxv)  = mxv;
	A(A < mnv)  = mnv;
end;

if ~isfield(V,'private') | ~isfield(V.private,'fid') | isempty(V.private.fid),
	mach = 'native';
	if swapped(dt),
		if spm_platform('bigend'), mach = 'ieee-le'; else, mach = 'ieee-be'; end;
	end; 
	[pth,nam,ext] = fileparts(V.fname);
	fname         = fullfile(pth,[nam, '.img']);
	fid           = fopen(fname,'r+',mach);
	if fid == -1,
		fid   = fopen(fname,'w',mach);
		if fid == -1,
			error(['Error opening ' fname '. Check that you have write permission.']);
		end;
	end;
else,
	if isempty(fopen(V.private.fid)),
		mach = 'native';
		if swapped(dt),
			if spm_platform('bigend'), mach = 'ieee-le'; else, mach = 'ieee-be'; end;
		end;
		V.private.fid = fopen(fname,'r+',mach);
		if V.private.fid == -1,
			error(['Error opening ' fname '. Check that you have write permission.']);
		end;
	end;
	fid = V.private.fid;
end;

% Seek to the appropriate offset
datasize = bits(dt)/8;
off   = (p-1)*datasize*prod(V.dim(1:2)) + V.pinfo(3,1);
fseek(fid,0,'bof'); % A bug in Matlab 6.5 means that a rewind is needed
if fseek(fid,off,'bof')==-1,
	% Need this because fseek in Matlab does not seek past the EOF
	fseek(fid,0,'bof'); % A bug in Matlab 6.5 means that a rewind is needed
	fseek(fid,0,'eof');
	curr_off = ftell(fid);
	blanks   = zeros(off-curr_off,1);
	if fwrite(fid,blanks,'uchar') ~= prod(size(blanks)),
		write_error_message(V.fname);
		error(['Error writing ' V.fname '.']);
	end;
	fseek(fid,0,'bof'); % A bug in Matlab 6.5 means that a rewind is needed
	if fseek(fid,off,'bof') == -1,
		write_error_message(V.fname);
		error(['Error writing ' V.fname '.']);
		return;
	end;
end;

if fwrite(fid,A,deblank(prec(dt,:))) ~= prod(size(A)),
	write_error_message(V.fname);
	error(['Error writing ' V.fname '.']);
end;

if ~isfield(V,'private') | ~isfield(V.private,'fid') | isempty(V.private.fid), fclose(fid); end;

return;
%_______________________________________________________________________

%_______________________________________________________________________
function write_error_message(q)
str = {...
	'Error writing:',...
	' ',...
	['        ',spm_str_manip(q,'k40d')],...
	' ',...
	'Check disk space / disk quota.'};
spm('alert*',str,mfilename,sqrt(-1));

return;
%_______________________________________________________________________
