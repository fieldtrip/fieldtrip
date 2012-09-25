function V = spm_vol_ana(fname, n)
% Get header information etc for Analyze 7.5 image
%
% FORMAT V = spm_vol_ana(fname, n)
% fname - name of Analyze file
% n     - frame number (defaults to 1)
% V     -  a vector of structures containing image volume information.
%          See spm_vol.m for more information
%
% Note that spm_flip_analyze_images.m should be configured for your
% data.  Images that were flipped at spatial normalisation in SPM99
% should be set so that they are left-right flipped.  Images that
% were not flipped at spatial normalisation should be set to be
% unflipped.
%
% Most fields in an Analyze .hdr are ignored, except...
%
% 	* hdr.dime.datatype is used to determine the datatype of the file.
% 	* hdr.dime.dim(1)<0 | hdr.dime.dim(1)>15 indicate that values in the
% 	  .hdr and .img should be flipped.
% 	* hdr.dime.dim(2) - hdr.dime.dim(5) are the x,y,z and t dimensions.
% 	* If the .hdr is not customised for storing orientation information
% 	  then the following fields may be read to derive the information:
% 		* .mat file if it exists
% 		* otherwise hdr.hist.origin to get origin of volume
% 		  and hdr.dime.pixdim(2:4) to get the voxel sizes.
%               * if hdr.hist.origin is set to [0 0 0], then the origin is
%                 assumed to be at the centre of the volume
%                 - (hdr.dime.dim(2:4)+1)/2.
% 	* Scalefactors and dc-offset are derived from hdr.dime.funused1
%         hdr.dime.funused2 respectively (if funused1~=0 & isfinite(funused1)).
% 	  If hdr.dime.funused1 is zero or non-finite then they are derived
%         from hdr.dime.cal_max, hdr.dime.cal_min, hdr.dime.glmax and
%         hdr.dime.glmin.
% 	* hdr.dime.vox_offset to get the offset into the volume.
%
%_______________________________________________________________________
% @(#)spm_vol_ana.m	2.7 John Ashburner (from Karl's old code) 03/05/27

if nargin<2, n = 1; end;
if ischar(n), n = str2num(n); end;

%
%-----------------------------------------------------------------------
[pth,nam,ext] = fileparts(deblank(fname));
hfname = fullfile(pth,[nam '.hdr']);
mfname = fullfile(pth,[nam '.mat']);
fname  = fullfile(pth,[nam   ext]);

if ~strcmp(ext,'.img'), error(['"' fname '" is not a ".img" file.']); end;

% Read the .hdr information
%-----------------------------------------------------------------------
[hdr,otherendian] = spm_read_hdr(hfname);
if isempty(hdr),
	V = [];
	return;
end;

% A little fix that sometimes helps
%-----------------------------------------------------------------------
if (hdr.dime.dim(1) < 3) & (hdr.dime.dim(4) == 0), hdr.dime.dim(4) = 1; end;
if (hdr.dime.dim(1) < 4) & (hdr.dime.dim(5) == 0), hdr.dime.dim(5) = 1; end;

if hdr.dime.dim(5)<n,
	error(['Not enough volumes in "' hfname '" (' num2str(n) '>' num2str(hdr.dime.dim(5)) ').']);
end;

% Datatype - set image to be byte-swapped if necessary
%-----------------------------------------------------------------------
dt  = hdr.dime.datatype;
if otherendian & dt~=2 & dt~=2+128, dt = dt*256; end;
dim     = [hdr.dime.dim(2:4) dt];
%-----------------------------------------------------------------------

% Orientation information...
%-----------------------------------------------------------------------
if hdr.hist.orient, warning(['Ignoring hdr.hist.orient field of "' hfname '".']); end;

mat = [];
if isfield(hdr,'spmf'),
	mat = hdr.spmf(n).mat;
else,
	if exist(mfname)==2,
		t   = load(mfname);
		if isfield(t,'mat'),
			if size(t.mat,3)>=n  & any(any(t.mat(:,:,n))),
				mat = t.mat(:,:,n);
			end;
		elseif isfield(t,'M'),
			mat = t.M;
			if spm_flip_analyze_images, mat = diag([-1 1 1 1])*mat; end;
		end;
	end;
	if isempty(mat),
		if any(hdr.hist.origin(1:3)),
			origin = hdr.hist.origin(1:3);
		else,
			origin = (hdr.dime.dim(2:4)+1)/2;
		end;
		vox    = hdr.dime.pixdim(2:4);
		if all(vox == 0), vox = [1 1 1]; end;
		off    = -vox.*origin;
		mat    = [vox(1) 0 0 off(1) ; 0 vox(2) 0 off(2) ; 0 0 vox(3) off(3) ; 0 0 0 1];
		if spm_flip_analyze_images, mat = diag([-1 1 1 1])*mat; end;
	end;
end;
%-----------------------------------------------------------------------

% Scaling etc... and offset into file...
%-----------------------------------------------------------------------
if isfinite(hdr.dime.funused1) & hdr.dime.funused1,
	scal  = hdr.dime.funused1;
	if isfinite(hdr.dime.funused2),
		dcoff = hdr.dime.funused2;
	else,
		dcoff = 0;
	end;
else,
	if hdr.dime.glmax-hdr.dime.glmin & hdr.dime.cal_max-hdr.dime.cal_min,
		scal  = (hdr.dime.cal_max-hdr.dime.cal_min)/(hdr.dime.glmax-hdr.dime.glmin);
		dcoff = hdr.dime.cal_min - scal*hdr.dime.glmin;
		%if hdr.dime.funused1,
		%	warning(['Taking scalefactor etc from glmax/glmin & cal_max/cal_min of "' fname '".']);
		%end;
	else,
		scal  = 1;
		dcoff = 0;
		warning(['Assuming a scalefactor of 1 for "' fname '".']);
	end;
end;

offset = hdr.dime.vox_offset + (n-1)*prod(dim(1:3))*spm_type(hdr.dime.datatype,'bits')/8;

if rem(offset, spm_type(hdr.dime.datatype,'bits')/8),
	error(['Offset into file must be a multiple of ' spm_type(hdr.dime.datatype,'bits')/8 '.']);
end;

pinfo   = [scal dcoff offset]';
%-----------------------------------------------------------------------

private = struct('hdr',hdr);

% Make volume handle
%-----------------------------------------------------------------------
descrip = hdr.hist.descrip;
V       = struct(...
	'fname',	fname,...
	'dim',		dim,...
	'pinfo',	pinfo,...
	'mat',		mat,...
	'descrip',	descrip,...
	'n',		n,...
	'private',	private);
%-----------------------------------------------------------------------
return;
%_______________________________________________________________________
%_______________________________________________________________________
