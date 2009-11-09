function V=spm_vol_minc(fname,n)
% Get header information etc. for MINC images.
%  FORMAT V = spm_vol(P)
%  P - a MINC filename.
%  V - a structure containing image volume information.
%
% The elements of V are described by the help for spm_vol, except for
% an additional field (V.private.cdf) that contains the NetCDF header
% information.
%
% The MINC file format was developed by Peter Neelin at the Montréal
% Neurological Institute, and is based upon the NetCDF libraries.
% The NetCDF documentation specifically recommends that people do not
% write their own libraries for accessing the data.  However, this
% was necessary in this case, in order to attempt to get MINC images into
% a format that the spm image reading utilities could use.
% _______________________________________________________________________
%  @(#)spm_vol_minc.m	2.16 John Ashburner 03/04/11

if nargin<2, n = 1; end;
if ischar(n), n = str2num(n); end;

cdf = spm_read_netcdf(fname);
if isempty(cdf), V=[]; return; end;

d_types     = [2 2 4 8 16 64];
dsizes      = [1 1 2 4  4  8];
space_names = {'xspace','yspace','zspace'};
img         = findvar(cdf.var_array,'image');
nd          = length(img.dimid);

if nd<3, error(['Not enough dimensions in "' fname '"']); end;
%for i=1:3,
%	if ~strcmp(space_names{i},cdf.dim_array(img.dimid(nd+1-i)).name),
%		error(['Incompatible dimension names in "' fname '"']);
%	end;
%end;

dim = zeros(1,nd);
for i=1:nd,
	dim(i) = cdf.dim_array(img.dimid(nd+1-i)).dim_length;
end;
dim = [dim 1];
if prod(dim(4:end))<n, error(['No volume # ' num2str(n) ' in "' fname '".']);end;

% Extract the dimensions.
%-----------------------------------------------------------------------
dim = dim(1:min(size(dim,2),3));
dim = [dim ones(1,3-size(dim,2))];

datatype = d_types(img.nc_type);
signed   = findvar(img.vatt_array,'signtype');
signed   = strcmp(signed.val,'signed__');
is_flt   = 0;
switch datatype,
	case {2},
		if signed,
			datatype = datatype + 128;
			range = [-2^7 2^7-1]';
		else,
			range = [0 2^8-1]';
		end;
	case {4},
		if signed,
			range = [-2^15 2^15-1]';
		else,
			datatype = datatype + 128;
			range = [0 2^16-1]';
		end;
	case {8},
		if signed,
			range = [-2^31 2^31-1]';
		else,
			datatype = datatype + 128;
			range = [0 2^32-1]';
		end;
	otherwise,
		is_flt = 1;
end;
if ~spm_platform('bigend') & datatype~=2 & datatype~=2+128, datatype = datatype*256; end;

dim   = [dim(1:3) datatype];

if ~is_flt,
	tmp = findvar(img.vatt_array,'valid_range');
	if isempty(tmp),
		disp(['Can''t get valid_range for "' fname '" - having to guess']);
	else
		range = tmp.val;
	end;

	fp   = fopen(fname,'r','ieee-be');
	imax = get_imax(fp, cdf, 'image-max', fname, 1, n, dim);
	imin = get_imax(fp, cdf, 'image-min', fname, 0, n, dim);
	fclose(fp);

	scale = (imax-imin)/(range(2)-range(1));
	dcoff = imin-range(1)*scale;
else,
	scale =  ones(1,dim(3));
	dcoff = zeros(1,dim(3));
end;

off   = img.begin + (n-1)*dsizes(img.nc_type)*prod(dim(1:3));
psize = dsizes(img.nc_type)*prod(dim(1:2));
off   = cumsum(ones(1,dim(3))*psize)-psize+off;

pinfo = [scale ; dcoff ; off];


% Extract affine transformation from voxel to world co-ordinates
%-----------------------------------------------------------------------
step  = [1 1 1];
start = [0 0 0]';
dircos = eye(3);

for j=1:3,
	nam    = cdf.dim_array(img.dimid(nd+1-j)).name;
	space  = findvar(cdf.var_array,nam);
	tmp    = findvar(space.vatt_array,'step');
	if ~isempty(tmp), step(j) = tmp.val; end;
	tmp    = findvar(space.vatt_array,'start');
	if ~isempty(tmp), start(j) = tmp.val; else, start(j) = -dim(j)/2*step(j); end;
	tmp    = findvar(space.vatt_array,'direction_cosines');
	if ~isempty(tmp), dircos(:,j) = tmp.val; end;
end;

shiftm = [1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1];
mat    = [[dircos*diag(step) dircos*start] ; [0 0 0 1]] * shiftm;

% Because there are not yet any routines to write the matrix information
% to MINC files, any changes to the matrix values will be made to `.mat'
% files.  The values in the `.mat' files should override the values from the
% headers.
matname = [spm_str_manip(fname,'sd') '.mat'];
if exist(matname) == 2,
	str=load(matname);
	if isfield(str,'mat'),
		if size(str.mat,3)>=n  & any(any(str.mat(:,:,n))),
			mat = str.mat(:,:,n);
		end;
	elseif isfield(str,'M'),
		mat = str.M;
	end;
end;

private = struct('cdf',cdf);
V       = struct('fname',fname,'dim',dim,'mat',mat,'pinfo',pinfo,'n',1,...
		'descrip','MINC file','private',private);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function var = findvar(varlist, name)
% Finds the structure in a list of structures that has a name element
% matching the second argument.
for i=1:prod(size(varlist)),
	if strcmp(varlist(i).name,name),
		var = varlist(i);
		return;
	end;
end;
var = [];
%error(['Can''t find "' name '".']);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function str = dtypestr(i)
% Returns a string appropriate for reading or writing the CDF data-type.
types = str2mat('uint8','uint8','int16','int32','float','double');
str   = deblank(types(i,:));
return;
%_______________________________________________________________________

%_______________________________________________________________________
function imax = get_imax(fp, cdf, strng, fname, def, n, dim)
img     = findvar(cdf.var_array,'image');
nd      = length(img.dimid);
alldims = zeros(1,length(cdf.dim_array));
for i=1:length(alldims),
	alldims(i) = cdf.dim_array(i).dim_length;
end;

str  = findvar(cdf.var_array,strng);

if ~isempty(str) & str.nc_type == 6,
	if any(str.dimid == img.dimid(nd)) | any(str.dimid == img.dimid(nd-1)),
		error(['Fastest two dimensions of "' fname '" should all have the same scalefactors etc']);
	end;

	ddim = fliplr(alldims(str.dimid));
	fseek(fp,str.begin,'bof');
	nel  = str.vsize/(spm_type(dtypestr(str.nc_type),'bits')/8);
	imax = fread(fp,nel,dtypestr(str.nc_type))';

	if nel==1,
		imax = imax*ones(1,dim(3));
	elseif nel== prod(ddim),
		imax = reshape(imax,[ddim 1 1]);

		% This bit is a fudge and may not work in all situations
		% ======================================================
		if str.dimid(end)==img.dimid(end-2),
			imax = imax(:,n)';
		else,
			imax = imax(n)*ones(1,dim(3));
		end;
		% ======================================================

	else,
		error(['Problem with ' strng]);
	end;
else,
	imax = ones(1,dim(3))*def;
	disp(['Can''t get ' strng ' for "' fname '" - guessing it is' num2str(def) '.']);
end;
return;

